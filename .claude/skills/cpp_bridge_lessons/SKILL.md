---
name: cpp_bridge_lessons
version: "0.4.0"
description: Reference material for wrapping a container-based MOM6 Fortran subroutine in a runtime-dispatched C++/AMReX bridge -- the shim and interface templates, Fortran/C type mapping, capture-record naming, module imports, and recurring pitfalls, in numbered sections §1-§17. Companion to generate_cpp_bridge, which cites these sections by number; invoke once per session before running it.
user-invocable: true
---

# C++/AMReX bridge reference

Invoke once per session before `generate_cpp_bridge`, which is the procedure and cites these
sections as "lessons.md §N". Sections are numbered §1–§17; never renumber them. The pattern comes
from TURBO-ESM/MOM6 PR #15 (`PPM_limit_pos`, `PPM_limit_cw84`, `PPM_reconstruction_y` in
`MOM_continuity_PPM`).

**Standing rule:** never reword, shorten or "fix" an existing Fortran comment; doc comments are
copied verbatim (`array_container_lessons`).

## 1. The pattern

For a container-based subroutine `FOO`:
1. The implementation is renamed `FOO_fortran`, otherwise byte-for-byte unchanged (§12).
2. A new `FOO` with the identical public signature becomes a dispatcher ("shim") choosing, per
   kernel at run time, Fortran, capture, or the C++/AMReX bridge (§2, §7).
3. A `bind(C)` interface declares `<prefix>_FOO_bridge`, implemented on the AMReX side (§3.3).
4. Containers (`RealArray_t`, `Box_t`) cross the boundary as POD mirrors (`RealArray_C`, `Box_C`)
   without copying data (§3.1).

Callers don't change. The default mode is Fortran, so a build with no env var set is bit-identical
to before.

## 2. The dispatcher shim

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

- The shim's signature is the converted subroutine's, verbatim, so callers don't change.
- One dedicated env var per kernel (§7); `case default` is the Fortran fallback for unknown modes
  and for builds without `_TIM`.
- Capture records inputs before and outputs after the Fortran truth (§6, §13), once per kernel per
  run, on the root PE only.

## 3. The C interoperability layer

### 3.1 Type mapping

| Fortran side | `bind(C)` side | Holds |
|---|---|---|
| `RealArray_t` | `RealArray_C` | n-D `real(real64)` array + shape/lb/ub/rank |
| `IntArray_t`  | `IntArray_C`  | n-D integer array + shape/lb/ub/rank |
| `Box_t`       | `Box_C`       | iteration range: two `c_ptr`s, to `idxS` and `idxE` (§5) |

```fortran
type, bind(C) :: RealArray_C
  type(c_ptr)    :: data        ! &this%data(1)
  type(c_ptr)    :: shape       ! &this%shape(1)
  type(c_ptr)    :: lb          ! &this%lb(1)
  type(c_ptr)    :: ub          ! &this%ub(1)
  integer(c_int) :: rank
end type
```

`%to_c()` fills the mirror. In TIM's infra `to_c_Real`/`to_c_Int` are null-safe, like
`Box_t%to_c`: an unassociated container (e.g. an absent optional) becomes an all-null descriptor
with `rank=0` instead of crashing on `c_loc(this%data(1))`, so the C side does the null check:

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

An older `array_mod.F90` may lack the null-safe form; check the checkout before relying on it, and
otherwise never call `%to_c()` on an unallocated container. FMS2's `to_c` is commented out
entirely, which is why every `%to_c()` must sit inside `#ifdef _TIM` (§8).

FMS2's `array_mod`/`box_mod` define the `_C` mirrors; TIM's import `RealArray_C`, `IntArray_C`,
`LogicalArray_C` and `Box_C` from TIM's `turbotmp_bridge_c_types` (C header
`turbotmp_bridge_c_types.h`) and re-export them, so `use array_mod, only : RealArray_C` works
under either. `bind(C)` matches components by position, not name (Fortran `rank` is C `dim`), so
TIM's module and header must agree in order and kind; nothing checks this at build or link time.

### 3.2 Scalars

Every scalar in a `bind(C)` interface is passed by value with an `iso_c_binding` kind:

| Fortran | `bind(C)` dummy |
|---|---|
| `real` | `real(c_double), intent(in), value` |
| `integer` | `integer(c_int), intent(in), value` |
| `logical` | `logical(c_bool), intent(in), value` (default `logical` is not `c_bool`) |
| pointer (e.g. `OBC`) | `type(c_ptr), intent(in), value` |

A default-kind `logical` passed for `c_bool` is a silent ABI bug, so copy it into a
`logical(c_bool)` local first; pass `int(n, c_int)` where default integer may be 8 bytes. Never
`c_loc` a disassociated pointer:

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

If the pointer dummy is also `optional` (MOM6's `OBC` often is), `associated()` on an absent one is
illegal: set a logical once at the shim's top, `OBC_associated = .false. ; if (present(OBC))
OBC_associated = associated(OBC)`, and use it for both the `c_loc` test and the capture record
`_OBC_associated`.

### 3.3 The interface block

One block per kernel at the top of the host module, every dummy doc-commented:

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

Name: `<prefix>_<original_name>_bridge`. Use the prefix the tree's existing `_bridge) bind(C)`
declarations use; otherwise `turbotmp_`, which marks the symbol as temporary scaffolding.

## 4. Call sites (background)

Wrapping raw arrays into containers at call sites is `convert_array_containers`'s job
(`array_container_lessons` is the API reference); the bridge never adds it. The recipe is box →
`%alloc(lb=…, ub=…, source=x)` (omit `source` for pure outputs) → call → `%copy2F` written ones →
`%free`. `%dup` no longer exists; `%copy2Array` only refills an already-allocated container.

## 5. `Box_t`

Kernels take a `Box_t` instead of loop-bound integers and loop with
`do concurrent (k=bx%idxS(3):bx%idxE(3), j=bx%idxS(2):bx%idxE(2), i=bx%idxS(1):bx%idxE(1))`.
Introducing it is `convert_array_containers`'s job (`array_container_lessons` §5). Its mirror
`Box_C` carries `idxS`/`idxE`, which AMReX consumes as a `Box`.

## 6. Capture mode

`turbotmp_helperF::io_recorder` writes paired files: `*.bin` (stream unformatted) and `*.meta`
(one line per entry: `name type byte-offset`). `add` is generic over `RealArray_t`, `Box_t`, real,
integer and logical; `get` reads back. A C++ test driver replays the inputs through the AMReX
kernel and diffs against the `_after` records (§13). The `already_recorded`/`mark_recorded`
registry is required, not optional: a kernel runs many times per step, and each records once per
process.

## 7. Mode selection

`mode = getenv_mode("<KERNEL>_MODE", default=TIMH_runFORTRAN)` uppercases the env var's value:

| Value | Constant (`turbotmp_helperF`) | Behaviour |
|---|---|---|
| `FORTRAN` | `TIMH_runFORTRAN` (103) | native Fortran (default) |
| `CAPTURE` | `TIMH_capture` (102) | record inputs+outputs around the Fortran truth |
| `AMREX` | `TIMH_runAMREX` (101) | call the C++/AMReX bridge (requires `_TIM`) |

Unknown values are FATAL. Per-kernel env vars let one kernel run under AMReX while the rest stay
Fortran, for bring-up and bisection; never a global mode variable.

## 8. `#ifdef _TIM`

Only the `case (TIMH_runAMREX)` arm is gated, selector and body together (gating only the body
leaves a reachable empty arm that silently does nothing). Capture is never gated: it needs only
the recorder and the Fortran truth. Without `_TIM` no C++ symbol is referenced and the link needs
no AMReX.

## 9. Files involved

| File | Role |
|---|---|
| `src/framework/turbotmp_helperF.F90` | mode parser, kernel registry, `io_recorder` |
| `config_src/infra/{FMS2,TIM}/box_mod.F90` | `Box_t` / `Box_C` |
| `config_src/infra/{FMS2,TIM}/array_mod.F90` | `RealArray_t`, `IntArray_t`, their `_C` mirrors, `to_c`, `view`, `alloc`, `copy2F`, `copy2Array`, binary I/O |
| the host module (e.g. `src/core/MOM_continuity_PPM.F90`) | `*_fortran` rename, shim, interface block |

A build selects exactly one infra layer, so FMS2 and TIM must offer the same public symbols
(TIM's re-exports the `_C` mirrors rather than defining them, §3.1). The C++ side (kernels,
`extern "C"` bridges, `turbotmp_helper`) is in protoMOMxx `src/mom_kernels/`, built as
`protoMOMxx::mom_bridge`; under `MOM6_INFRA=TIM`, `src/CMakeLists.txt` finds and links it.

## 10. Recurring pitfalls

1. **Kinds at the boundary.** Default `logical`/`integer` passed for `c_bool`/`c_int` (§3.2).
2. **Null actuals.** `c_loc` on a disassociated pointer, or `%to_c()` on an unallocated container
   where `to_c` is not null-safe (§3.1).
3. **Doc comments.** Every dummy of the shim and the interface block needs `!<` or `!!`; write
   them first.
4. **`#ifdef _TIM` around the body only**, not the `case` selector (§8).
5. **Incomplete flattening of a derived-type dummy.** A struct crosses `bind(C)` as one channel per
   field, and a hand-built field list is easy to get short. Derive it from this subroutine's own
   `_fortran` body by grep, including call actuals and callees the struct is passed to, never from
   a sibling bridge on the same type; give a field behind its own guard its own null-able channel;
   check the field sets agree by script. Nothing else catches it: FORTRAN mode is unchanged,
   argument counts agree, and the capture list omits the same field. Example:
   `zonal_mass_flux`'s bridge copied `set_zonal_BT_cont`'s six fields and missed `BT_cont%h_u`
   (written once, as a call actual, behind `%associated()`), leaving it stale in AMREX mode.
6. **Name clashes are case-insensitive.** A local `coupling_state_c` collides with a module type
   `coupling_state_C` (gfortran: "cannot be host associated"). Check every new local and dummy
   against the module's type and procedure names, not only the routine's; e.g. use
   `<name>_cdesc`.
7. **Capture without the registry** writes one record per invocation (§6).

## 11. Invariants

- The public signature of `FOO` (callers don't change).
- `FOO_fortran` is the original, renamed only; default-mode numerics are bit-identical.
- The Fortran and AMReX kernels consume the same `Box_t` range and the same array bounds; capture
  mode is what validates this.

## 12. Renaming the original to `*_fortran`

The rename of the `subroutine`/`end subroutine` lines is the only edit. A correctly converted
subroutine already looks like this; if it doesn't (raw array dummies, loop-bound integers instead
of a `Box_t`), stop and run `convert_array_containers` first rather than converting inline.

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

## 13. Capture-record naming

| Argument | Record name |
|---|---|
| read only | `_<argname>` |
| read and written | `_<argname>_before`, then `_<argname>_after` |
| written only | `_<argname>_after` |
| iteration box | `_<bx>` |
| scalar | `_<argname>` (no before/after) |

The `kernel` string is short, lowercased and unique per kernel; it names
`capture/<kernel>.bin` and `capture/<kernel>.meta`.

## 14. Module imports for the shim

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

Deduplicate against the module's existing `use`s; `is_root_pe` usually already comes from
`MOM_error_handler`. Never re-implement anything `turbotmp_helperF` provides.

## 15. Halo checks belong outside the kernel

A halo-sufficiency `MOM_error(FATAL, ...)` inside the kernel moves to the caller or the shim
entry, so it runs in every mode and the C++ side can assume sufficient halos.

## 16. CPU clocks stay at the caller

`cpu_clock_begin`/`cpu_clock_end` stay around the shim call at the caller, never inside the shim
(that would double-clock and hide per-mode timing).

## 17. Verification matrix

| Mode | Setting | Expected |
|---|---|---|
| FORTRAN | env var unset | bit-identical to the pre-bridge build |
| CAPTURE | `<ENV_VAR>=CAPTURE` | `capture/<kernel>.{bin,meta}` exist, one record per kernel, from rank 0 only |
| AMREX | `<ENV_VAR>=AMREX`, built with `_TIM` | replay diff against the CAPTURE outputs is zero (or within a stated tolerance) |

Building and running these is the user's job unless they explicitly allow it in this
conversation. If only the shim and capture are delivered, the C++ side of
`<prefix>_<kernel>_bridge` is the next deliverable.
