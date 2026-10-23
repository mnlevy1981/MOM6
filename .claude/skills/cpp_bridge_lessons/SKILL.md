---
name: cpp_bridge_lessons
version: "0.3.2"
description: Reference material for wrapping a MOM6 Fortran subroutine in a runtime-dispatched C++/AMReX bridge -- type-mapping tables, code templates, naming conventions, and recurring pitfalls, organized in numbered sections §1-§17. Companion to the generate_cpp_bridge skill, which cites these sections by number (e.g. "lessons.md §12") throughout its procedure. Invoke once near the start of a session that will run generate_cpp_bridge one or more times; not needed again per-subroutine within that session.
user-invocable: true
---

# Lessons from PR #15: AMReX/C++ Bridge for PPM Subroutines

> Invoke once, near the start of a session, before running
> `generate_cpp_bridge` (or `generate_amrex_code`). Stays available for
> the rest of the session. Sections are numbered §1–§17 and cited by
> that numbering elsewhere; do not renumber them.

Source: [TURBO-ESM/MOM6 PR #15](https://github.com/TURBO-ESM/MOM6/pull/15)
— merged commit `daf6abefb`. This document distills the design and
patterns it used to wrap three subroutines in `MOM_continuity_PPM`
(`PPM_limit_pos`, `PPM_limit_cw84`, `PPM_reconstruction_y`) so they can
be redirected at runtime to C++/AMReX, while keeping the Fortran truth
available and adding a capture mode for offline validation. It's a
*non-invasive* refactor — callers don't change; the subroutine they call
becomes a small dispatcher (a "shim") picking between three implementations.

**Pre-condition:** work in a TURBO-ESM/MOM6 checkout with `src/` and
`config_src/`, on (or rebased onto) `dev/turbo-debug`. Cloning isn't
performed —
`git clone -b dev/turbo-debug git@github.com:TURBO-ESM/MOM6.git <dir>`
once, then pass that directory as `<work-directory>` on every invocation.

---

## 1. The big picture — what the PR actually does

For each target subroutine `FOO`:

1. The original implementation is **renamed** `FOO_fortran` and left in
   place, byte-for-byte equivalent to before.
2. A new subroutine **`FOO`** is introduced with the *same public signature*.
   It is a thin dispatcher that selects one of three execution modes via an
   environment variable.
3. A `bind(C)` interface block declares an opaque **`<prefix>_FOO_bridge`**
   that is implemented in C++ on the AMReX side.
4. New container types (`RealArray_t`, `Box_t`) and their `bind(C)` mirrors
   (`RealArray_C`, `Box_C`) carry array data and iteration ranges across the
   Fortran↔C++ boundary without any data copy.
5. A small framework module ([`turbotmp_helperF.F90`](../../../src/framework/turbotmp_helperF.F90))
   provides the mode parser, the kernel registry, and a binary/metadata
   `io_recorder` for capture mode.
6. Every existing caller of `FOO` is rewritten to wrap its raw Fortran arrays
   in `RealArray_t` containers, call the shim, and copy results back.

Because the shim's default mode is `TIMH_runFORTRAN`, builds with no env var
set produce *bit-identical* numerics to the pre-PR code. The CESM regression
test confirmed zero numerical change.

---

## 2. The dispatcher pattern (the "shim")

```fortran
!< shim for PPM_limit_pos
subroutine PPM_limit_pos(bx, h_in, h_L, h_R, h_min)
  type(box_t),       intent(in)    :: bx
  type(RealArray_t), intent(in)    :: h_in
  type(RealArray_t), intent(inout) :: h_L
  type(RealArray_t), intent(inout) :: h_R
  real,              intent(in)    :: h_min

  integer            :: mode
  type(RealArray_C)  :: h_in_c, h_L_c, h_R_c
  type(Box_c)        :: bx_c
  type(io_recorder)  :: rec
  logical            :: capture
  character(len=80)  :: kernel
  character(len=256) :: binFile, metaFile

  kernel = "ppm_limit_pos"
  mode = getenv_mode("PPM_LIMIT_POS_MODE", default=TIMH_runFORTRAN)

  select case (mode)
    case (TIMH_capture)
       capture = (.not. already_recorded(trim(kernel))) .and. is_root_pe()

       if (capture) then
         call rec%open_write("capture/"//trim(kernel)//".bin", &
                             "capture/"//trim(kernel)//".meta")
         call rec%add("_bx",         bx)
         call rec%add("_h_in",       h_in)
         call rec%add("_h_L_before", h_L)
         call rec%add("_h_R_before", h_R)
         call rec%add("_h_min",      h_min)
       endif

       call ppm_limit_pos_fortran(bx, h_in, h_L, h_R, h_min)   ! truth

       if (capture) then
         call rec%add("_h_L_after", h_L)
         call rec%add("_h_R_after", h_R)
         call rec%close()
         call mark_recorded(trim(kernel))
       endif

#ifdef _TIM
    case (TIMH_runAMREX)
       bx_c   = bx%to_c()
       h_in_c = h_in%to_c() ; h_L_c = h_L%to_c() ; h_R_c = h_R%to_c()
       call turbotmp_ppm_limit_pos_bridge(bx_c, h_in_c, h_L_c, h_R_c, h_min)
#endif

    case default
       call ppm_limit_pos_fortran(bx, h_in, h_L, h_R, h_min)
  end select
end subroutine PPM_limit_pos
```

**Key properties:**

- The shim's *own* signature is **identical** to the public signature the
  caller already used, so callers do not change.
- `getenv_mode("PPM_LIMIT_POS_MODE", default=TIMH_runFORTRAN)` — each kernel
  has a **dedicated env var**, so kernels can be enabled independently.
- `case default` falls back to Fortran; this is the safety net for unknown
  modes and for builds without `_TIM` defined.
- The C++ branch is wrapped in `#ifdef _TIM` so the bridge symbols are not
  required when AMReX/C++ is not linked.
- `capture` mode runs the Fortran truth, recording inputs *before* and
  outputs *after*. The "before/after" naming is what lets a C++-side test
  driver replay inputs and diff outputs.
- Capture is gated by `already_recorded(kernel) .and. is_root_pe()` — one
  record per kernel per run, on the root rank only.

---

## 3. The C interoperability layer

### 3.1 Two parallel type families

| Fortran-side (rich) | C-side (POD, `bind(C)`) | Purpose |
|---|---|---|
| `RealArray_t`       | `RealArray_C`           | n-D real(real64) array + shape/lb/ub/rank |
| `IntArray_t`        | `IntArray_C`            | n-D integer array  + shape/lb/ub/rank     |
| `Box_t`             | `Box_C`                 | iteration index range (idxS, idxE)        |

The Fortran-side type carries pointers, type-bound procedures, and *also*
serves as the value type passed in capture mode. The C-side struct is a
fixed-layout POD of `c_ptr`s + an `integer(c_int) :: rank` so it can be
declared identically in C++.

Where the C-side types are defined depends on the infra layer:

- **FMS2:** `config_src/infra/FMS2/array_mod.F90` and `box_mod.F90` define
  them with `type, bind(C)` blocks like the one below. FMS2 builds have no
  TIM.
- **TIM:** `config_src/infra/TIM/array_mod.F90` and `box_mod.F90` define
  none of them. They `use turbotmp_bridge_c_types` from TIM
  (`tim/fortran/turbotmp_bridge_c_types.F90`), which defines `RealArray_C`,
  `IntArray_C`, `LogicalArray_C` and `Box_C`, and re-export them in their
  `public` lists. The C++ side includes the matching header,
  `tim/cpp/turbotmp_bridge_c_types.h`.

Both layers offer the same public names, so
`use array_mod, only : RealArray_C` works on either. TIM's Fortran module and
its C header must agree in component order and type for every struct:
`bind(C)` matches components by position, not by name (the Fortran `rank` is
the C `dim`), and nothing checks this at build or link time.

```fortran
type, bind(C) :: RealArray_C
  type(c_ptr)    :: data        ! &this%data(1)
  type(c_ptr)    :: shape       ! &this%shape(1)
  type(c_ptr)    :: lb          ! &this%lb(1)
  type(c_ptr)    :: ub          ! &this%ub(1)
  integer(c_int) :: rank
end type
```

The Fortran-side type has a `to_c()` method that fills a `RealArray_C`.
In TIM's infra, `to_c_Real`/`to_c_Int` are **null-safe**: an unassociated
`data` pointer (e.g. an absent optional container) converts to an
all-null `cdesc` with `rank=0`, rather than crashing on
`c_loc(this%data(1))`. This mirrors `Box_t%to_c`'s existing pattern:

```fortran
function to_c_Real(this) result(cdesc)
  class(RealArray_t), intent(in) :: this
  type(RealArray_C) :: cdesc

  cdesc%data  = c_null_ptr
  cdesc%shape = c_null_ptr
  cdesc%lb    = c_null_ptr
  cdesc%ub    = c_null_ptr
  cdesc%rank  = 0

  if (associated(this%data)) then
    cdesc%data  = c_loc(this%data(1))
    cdesc%shape = c_loc(this%shape(1))
    cdesc%lb    = c_loc(this%lb(1))
    cdesc%ub    = c_loc(this%ub(1))
    cdesc%rank  = this%rank
  endif
end function
```

Bridge code passing an optional array container through `%to_c()` can
rely on the C-side null check rather than pre-verifying `associated()`
before calling `%to_c()` itself. FMS2's `to_c` is fully commented out
(no type-bound binding at all) and doesn't need this fix — every
`%to_c()` call site in the tree is gated by `#ifdef _TIM`, so FMS2
builds never compile or reach this code path.

### 3.2 Scalar passing rules (gotcha-prone)

All scalars in `bind(C)` interfaces are passed **by value** with explicit
`iso_c_binding` kinds:

- `real`        → `real(c_double),  intent(in), value`
- `integer`     → `integer(c_int),  intent(in), value`
- `logical`     → `logical(c_bool), intent(in), value`  (note: NOT default `logical`)
- generic pointer → `type(c_ptr),    intent(in), value`

At the call site, default-kind Fortran scalars must be cast:

```fortran
logical(c_bool) :: monotonic_c, simple_2nd_c
type(c_ptr)     :: OBC_c

monotonic_c  = monotonic
simple_2nd_c = simple_2nd
if (associated(OBC)) then
  OBC_c = c_loc(OBC)
else
  OBC_c = c_null_ptr
endif
```

### 3.3 Bridge interface declaration

The `bind(C)` interface block lives at the top of the host module
(`MOM_continuity_PPM`), one block per kernel:

```fortran
interface
  subroutine turbotmp_ppm_limit_pos_bridge(bx, h_in, h_L, h_R, h_min) bind(C)
    use iso_c_binding
    use array_mod, only : RealArray_c
    use box_mod,   only : Box_c
    implicit none
    type(Box_C),       intent(in)        :: bx
    type(RealArray_C), intent(in)        :: h_in
    type(RealArray_C), intent(inout)     :: h_L
    type(RealArray_C), intent(inout)     :: h_R
    real(c_double),    intent(in), value :: h_min
  end subroutine turbotmp_ppm_limit_pos_bridge
end interface
```

Naming convention: **`<prefix>_<original_name>_bridge`**. The PR uses the
`turbotmp_` prefix to mark the bridge as TURBO-tmp scaffolding — this is
deliberate: the prefix advertises that these symbols are temporary and may
be retired once full ports are validated.

---

## 4. Caller-side wrapping pattern

> **Scope note.** Converting raw arrays to containers is
> `convert_array_containers`'s job, not this one — run it first; this
> section is background so the bridge work reads sensibly, and because a
> call site may still be unconverted when its callee is bridged.
> Authoritative container API: `array_container_lessons`.

Call sites that invoke a container-based routine on raw arrays wrap those
arrays into `RealArray_t` containers. The recipe:

```fortran
type(RealArray_t) :: h_in_a, h_W_a, h_E_a
type(Box_t)       :: bx

call bx%safe_alloc(ndims=3)
call bx%set(idxS=[ish,jsh,1], idxE=[ieh,jeh,nz])

! Allocate to the array's bounds and copy the data in, in one call.
call h_in_a%alloc(lb=LBOUND(h_in), ub=UBOUND(h_in), source=h_in)
! Pure outputs: allocate only -- omit `source` rather than reading
! undefined memory.
call h_W_a%alloc(lb=LBOUND(h_W), ub=UBOUND(h_W))
call h_E_a%alloc(lb=LBOUND(h_E), ub=UBOUND(h_E))

if (monotonic) then
  call PPM_limit_cw84(bx, h_in_a, h_W_a, h_E_a)
else
  call PPM_limit_pos (bx, h_in_a, h_W_a, h_E_a, h_min)
endif

! Copy results back to native Fortran arrays
call h_W_a%copy2F(h_W)
call h_E_a%copy2F(h_E)

! Cleanup
call h_in_a%free() ; call h_W_a%free() ; call h_E_a%free()
call bx%free()
```

The pattern is: **alloc box → alloc containers → shim call → copy2F →
free**.

`%dup` no longer exists. It was restructured into `%alloc`: what was once
`call a%dup(x) ; call a%copy2Array(x)` is now the single call
`call a%alloc(lb=LBOUND(x), ub=UBOUND(x), source=x)`. Code emitted
against the old two-call recipe will not compile. `%copy2Array` still
exists and is the right tool for *re-filling* an already-allocated
container, but it does not allocate.

---

## 5. Iteration domain as a first-class object (`Box_t`)

> **Scope note.** Introducing a `Box_t` is `convert_array_containers`'s
> job (`array_container_lessons` §5 is the authoritative API reference);
> by the time this skill runs the kernel already takes one. Background
> plus the C-side mirror, which *is* bridge-specific.

Kernels take a `Box_t` rather than packing loop bounds into many integer
args, so loops inside them read:

```fortran
do concurrent (k=bx%idxS(3):bx%idxE(3), &
               j=bx%idxS(2):bx%idxE(2), &
               i=bx%idxS(1):bx%idxE(1))
  ...
enddo
```

`Box_t` supports `grow(dim, n)` and `shrink(dim, n)` to derive halo-extended
or halo-restricted boxes — these replace ad-hoc index arithmetic at the
caller and make the iteration domain explicit and inspectable.

The C-side mirror `Box_C` is just two `c_ptr`s, one each for `idxS` and
`idxE`. AMReX consumes the same start/end vectors as a `Box`.

---

## 6. Capture mode as a portable validation harness

`turbotmp_helperF::io_recorder` provides a **dump-and-replay** facility.

- `open_write` / `open_read` create paired `*.bin` (raw binary, stream
  unformatted) and `*.meta` (formatted, one entry per line: `name type
  byte-offset`) files.
- `add` is generic over `RealArray_t`, `Box_t`, real, integer, logical;
  `get` is the symmetric reader.
- The kernel writes `_input1`, `_input2`, … *then* runs the Fortran truth,
  *then* writes `_output1_after`, `_output2_after`. The C++ test driver
  reads the inputs, calls the AMReX kernel, and diffs against the
  `_after` arrays.

`already_recorded` / `mark_recorded` keep a `character(len=128) ::
recorded(max_kernels)` registry so each kernel records exactly once per
process — important because a single timestep may invoke a kernel hundreds
of times.

---

## 7. Mode selection — environment-variable driven

```fortran
mode = getenv_mode("PPM_LIMIT_POS_MODE", default=TIMH_runFORTRAN)
```

`getenv_mode` reads the named env var, uppercases it, and maps:

| String    | Constant            | Behaviour                                 |
|-----------|---------------------|-------------------------------------------|
| `FORTRAN` | `TIMH_runFORTRAN`   | Run native Fortran (default)              |
| `CAPTURE` | `TIMH_capture`      | Record inputs+outputs, run Fortran truth  |
| `AMREX`   | `TIMH_runAMREX`     | Call C++/AMReX bridge (requires `_TIM`)   |

Unknown values are FATAL. The integer constants live in `turbotmp_helperF`
(`TIMH_runAMREX = 101`, `TIMH_capture = 102`, `TIMH_runFORTRAN = 103`).

This is intentionally per-kernel: you can run one kernel under AMReX while
the rest stay on Fortran, simplifying bring-up and bisection.

---

## 8. Build-time gate on AMReX path (`#ifdef _TIM`)

Only the `TIMH_runAMREX` arms are wrapped in `#ifdef _TIM`. Capture mode is
*not* gated — it has no AMReX dependency, only the recorder + the existing
Fortran truth. This means:

- Pure-Fortran builds: omit `_TIM`, the C++ symbols are never referenced,
  link succeeds with no AMReX present.
- AMReX-enabled builds: define `_TIM`, link the C++ bridge library.
- Capture mode is always available.

---

## 9. Module/file layout introduced or changed

| File | Role |
|---|---|
| `src/framework/turbotmp_helperF.F90` (new) | mode parser, kernel registry, `io_recorder` |
| `config_src/infra/FMS2/box_mod.F90` (new) | `Box_t` / `Box_C` |
| `config_src/infra/TIM/box_mod.F90`  (new) | TIM-infra mirror of `box_mod` |
| `config_src/infra/FMS2/array_mod.F90` (extended) | `RealArray_t`, `IntArray_t`, `RealArray_C`, `IntArray_C`, `to_c`, `view*`, `alloc*`, `copy2F*`, `copy2Array*`, `dup*`, `write_binary`, `read_binary` |
| `config_src/infra/TIM/array_mod.F90`  (extended) | TIM-infra mirror of `array_mod` |
| `src/core/MOM_continuity_PPM.F90` (modified) | `*_fortran` rename, three new shims, three `bind(C)` interface blocks, callers rewritten |

The TIM and FMS2 infras keep parallel implementations because a build
selects exactly one infra layer; both must offer the same public symbols.
TIM's `array_mod`/`box_mod` still offer `RealArray_C`, `LogicalArray_C` and
`Box_C`, but no longer define them. They import them from TIM's
`turbotmp_bridge_c_types` (§3.1). FMS2's define their own.

The C++ side (kernels, `extern "C"` bridges, `turbotmp_helper`) lives in
protoMOMxx under `src/mom_kernels/`, built as `protoMOMxx::mom_kernels` and
`protoMOMxx::mom_bridge`. Under `MOM6_INFRA=TIM`, `src/CMakeLists.txt` finds
it with `find_package(protoMOMxx REQUIRED COMPONENTS MOM_BRIDGE)` and links
`protoMOMxx::mom_bridge` into `mom6_ocean`.

---

## 10. Recurring pitfalls observed in the commit log

The PR's commit history (~70 commits) reveals the friction worth flagging:

1. **`logical` → `logical(c_bool)` mismatches.** Fortran default `logical`
   is not `c_bool`; passing one for the other is a silent ABI bug. Always
   declare a `logical(c_bool) :: x_c` local and assign before the call.
2. **Default integers vs `c_int`.** Some platforms have 8-byte default
   integer; pass `int(iis, c_int)` if there is any doubt.
3. **Doxygen on shims.** Every dummy argument in the shim needs a `!<` or
   `!!` doc comment, including the bridge interface blocks. The PR has many
   "doxstring" fix commits — write the docs first.
4. **Null pointer args.** Optional pointer args (e.g., `OBC`) need
   `c_null_ptr` when not associated; never call `c_loc` on a disassociated
   pointer.
5. **`to_c` on an unassociated container.** TIM's `to_c_Real`/`to_c_Int`
   are null-safe (see §3.1) — an unassociated `data` pointer converts to
   an all-null `cdesc` (`rank=0`) instead of crashing on
   `c_loc(this%data(1))`. Don't assume this fix is present in every
   checkout, though: an older/unpatched copy of `array_mod.F90` still
   crashes on an unallocated container, so verify before relying on it.
6. **Capture file proliferation.** Without `already_recorded`, capture mode
   would produce one file per call site invocation. The registry is
   essential, not optional.
7. **`#ifdef _TIM` placement.** Wrap the *case arm body and the `case (...)`
   selector together*, not just the body — otherwise the unused `case` is
   reachable and falls through to fortran with no diagnostics. The PR's
   final form keeps both inside `#ifdef`.
8. **Incomplete flattening of a derived-type dummy.** A struct like
   `BT_cont` crosses `bind(C)` as one `RealArray_C` per field, so the
   field list is hand-built — and easy to get short. Real case:
   `zonal_mass_flux`'s bridge carried the six `FA_u_*`/`uBT_*` fields but
   not `BT_cont%h_u`, which `zonal_mass_flux_fortran` writes by passing
   it to `zonal_flux_thickness` at the very end of the body. Why it was
   missed: the list matched `set_zonal_BT_cont`'s bridge (correct *there*
   — that body never touches `h_u`); the six fields are `%view`'d in one
   block near the top while `h_u` appears only once, ~280 lines later,
   as a call argument; and `h_u` sits behind a second guard
   (`BT_cont%h_u%associated()`, allocated only with `alloc_faces`).
   Nothing caught it: FORTRAN mode stays bit-identical, argument counts
   between shim and interface agree, and CAPTURE omitted `h_u` from the
   same list, so a C++ replay passes. In AMREX mode `h_u` was simply
   left stale — the barotropic solver's face-thickness weights, silently
   wrong. The same gap propagated to `continuity_PPM`'s bridge.
   Fix: derive the list from the body by grep (including call actuals),
   give nested-optional fields their own null-able channel, and diff
   field sets by script (`generate_cpp_bridge` Steps 2 and 9).
9. **A shim's `_c` local named like a mirror type.** Fortran ignores
   case, so naming the C copy of a dummy `<dummy>_c` breaks whenever the
   module also defines a `bind(C)` mirror called `<dummy>_C`. Real case:
   in `MOM_vert_friction`, the Wave 1 shims declared
   `type(coupling_state_C) :: coupling_state_c`, and gfortran rejected
   it ("The type 'coupling_state_c' cannot be host associated ...
   blocked by an incompatible object of the same name"). Fix: give such a
   local a different name (e.g. `coupling_state_cdesc`), and check every
   new local and dummy case-insensitively against the module's type and
   procedure names, not only against the routine's own declarations.

---

## 11. What stays invariant when applying this pattern

- Public signature of the shimmed subroutine (so callers do not change).
- The original Fortran implementation, renamed to `*_fortran` and otherwise
  untouched. *Numerical truth must remain bit-identical* in the default
  mode.
- The kernel's iteration semantics: the Fortran kernel and the AMReX kernel
  must consume the same `Box_t` index range with the same lb/ub on the
  arrays — capture mode validates this.

---

## 12. Renaming the original to `*_fortran`

The original implementation is renamed `FOO` → `FOO_fortran`. **The rename
is the only edit.** Everything else below — container dummies, the
`%view` pointers, the `do concurrent` over `bx%idxS`/`bx%idxE` — is
already present, because `convert_array_containers` put it there before
this skill ran. The listing is here so you can recognise a correctly
converted subroutine, not as a transformation to perform.

If the subroutine does *not* already look like this, it has not been
converted: stop and run `convert_array_containers` on it first, rather
than converting it inline (see `array_container_lessons` §2).

```fortran
subroutine FOO_fortran(bx, h_in_a, h_L_a, h_R_a, h_min)
  type(Box_t),       intent(in)    :: bx
  type(RealArray_t), intent(in)    :: h_in_a
  type(RealArray_t), intent(inout) :: h_L_a
  type(RealArray_t), intent(inout) :: h_R_a
  real,              intent(in)    :: h_min

  real, dimension(:,:,:), contiguous, pointer :: h_in, h_L, h_R
  integer :: i, j, k

  call h_in_a%view(h_in)
  call h_L_a %view(h_L)
  call h_R_a %view(h_R)

  do concurrent (k=bx%idxS(3):bx%idxE(3), &
                 j=bx%idxS(2):bx%idxE(2), &
                 i=bx%idxS(1):bx%idxE(1))
    ! ...original body, unchanged...
  enddo
end subroutine FOO_fortran
```

If the original took loop-bound integers (`is, ie, js, je, nz`), collapse
them into one `Box_t` parameter; the explicit loop becomes the
`do concurrent` form above.

---

## 13. Capture-record naming convention

The recorder uses positional-but-readable names:

- input args: `_<argname>` (or `_<argname>_before` if also written as output)
- output args after the call: `_<argname>_after`
- iteration box: `_<bx>`
- scalars: `_<argname>` (no before/after — they don't change)

Choose a short lowercased `kernel` string per kernel; it becomes the
filename for `capture/<kernel>.bin` and `capture/<kernel>.meta`. Do not
reuse a kernel string across kernels.

---

## 14. Module `use` statements expected by the shim

The host module's preamble must import the bridge plumbing. A minimal set:

```fortran
use array_mod,         only : RealArray_t, RealArray_C
use box_mod,           only : Box_t, Box_C
use iso_c_binding,     only : c_double, c_int, c_ptr, c_loc, c_bool, &
                              c_null_char, c_null_ptr
use posix,             only : mkdir_posix
use turbotmp_helperF,  only : getenv_mode, io_recorder, &
                              already_recorded, mark_recorded
use turbotmp_helperF,  only : TIMH_runAMREX, TIMH_capture, TIMH_runFORTRAN
```

`is_root_pe` is typically already imported via `MOM_error_handler` — verify
before adding a duplicate `use`.

---

## 15. Halo-sufficiency checks belong outside the kernel

Halo-sufficiency `MOM_error(FATAL,...)` checks that previously sat *inside*
the kernel must move to the caller (or stay at the start of the shim) so
they execute regardless of mode. PR #15 relocated the halo checks for
`PPM_reconstruction_y` to its caller — this is what allows the C++ branch
to assume the halos are already big enough.

---

## 16. CPU clocks stay around the call site, not inside the shim

If `FOO_fortran` was previously inside a `cpu_clock_begin` /
`cpu_clock_end` pair, leave the clock calls at the *caller*, around the
shim invocation. Putting them inside the shim would double-clock when the
caller is also clocked, and would obscure mode-specific timing.

---

## 17. Verification matrix

After wrapping, run all three modes and verify:

| Mode      | Env var setting             | Expected outcome                                          |
|-----------|-----------------------------|-----------------------------------------------------------|
| FORTRAN   | unset (default)             | Bit-identical numerics vs. pre-PR build. CESM regression. |
| CAPTURE   | `<ENV_VAR>=CAPTURE`         | `capture/<kernel>.{bin,meta}` exists, exactly one record per kernel, written by rank 0 only. |
| AMREX     | `<ENV_VAR>=AMREX`, `_TIM`   | Replay-diff vs. CAPTURE outputs is zero (or within stated tolerance). |

If only the Fortran shim and capture mode are being delivered (no C++ side
yet), stop after CAPTURE verification and report that the AMReX-side
implementation of `<prefix>_<kernel>_bridge` is the next deliverable.
