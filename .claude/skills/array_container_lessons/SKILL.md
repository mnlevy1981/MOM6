---
name: array_container_lessons
version: "0.3.2"
description: Reference material for converting MOM6 subroutines from raw Fortran arrays (real, dimension(SZI_(G),SZJ_(G),SZK_(GV))) to the RealArray_t / IntArray_t / LogicalArray_t / Box_t container types -- the in-place conversion pattern, where the marshalling lives at the call site, the complete alloc/allocView/view/copy2F/copy2Array/free API with exact signatures, intent-driven copy rules, and recurring pitfalls, organized in numbered sections §1-§10. Companion to the convert_array_containers skill, which cites these sections by number throughout its procedure. Invoke once near the start of a session that will run convert_array_containers one or more times; not needed again per-subroutine within that session. STANDING RULE for every MOM6 edit: never reword, shorten or "fix" an existing Fortran comment; copy doc comments verbatim (see the box at the top).
user-invocable: true
---

# Array-container conversion: reference

> **STANDING RULE — Fortran comments are not yours to edit (user, repeated; 2026-09-29).**
> - Never change the wording of an existing Fortran comment: don't shorten it, paraphrase it,
>   drop an article ("the"), drop the trailing period, or "correct" a typo. This applies to `!<`,
>   `!!` and plain `!` comments alike.
> - When a new declaration takes its doc comment from somewhere else (e.g. a `grid_core` member
>   from `ocean_grid_type` in `MOM_grid.F90`, or a promoted dummy from the field it replaces),
>   **copy the text word for word**. Keep the upstream wording even where it looks wrong; for
>   example `IdyCv_OBCmask` says "1/dxCv".
> - **Line too long?** Move line breaks only: continue on `!!` lines aligned under the `!<`. Never
>   reword to make it fit. Keep a bracketed unit, and unit alternatives like `[degrees_N] or [km]
>   or [m]`, on one line.
> - **Verify before handing over:** extract the source comment by script and check word-for-word
>   equality with what you wrote.
> - **If a change makes an existing comment inaccurate** (e.g. a count or a rule it states), leave
>   the comment alone and tell the user.
> - Adding a brand-new comment for a brand-new declaration is fine. Editing any existing comment
>   is not.

> **STANDING RULE — use the pre-written skills (user, 2026-10-23).** Every MOM6 conversion step
> that a skill in `.claude/skills/` covers is done by invoking that skill with the Skill tool and
> following its procedure. Never do it by hand or with an ad-hoc editing script. A campaign is
> driven by invoking `convert_calltree` itself, not by reading its file. Hand-authored code is
> allowed only where `convert_calltree` Step 5a says so (Stage 1). See the standing rule at the
> top of `convert_calltree`.

> Invoke once, near the start of a session, before running
> `convert_array_containers`; it stays available for the rest of the
> session. Sections are numbered §1–§10 and cited by that number from
> the calling skill — do not renumber them.

Worked reference in the tree, `src/core/MOM_continuity_PPM.F90`:

- **A converted callee:** `zonal_edge_thickness_fortran` (~603) and its
  mirror `meridional_edge_thickness_fortran` (~747) — container dummies,
  `%view` pointers, `do concurrent` over a `Box_t`.
- **A call site doing the marshalling:** `zonal_edge_thickness` (~636)
  and `meridional_edge_thickness` (~780) — `alloc`/`copy2F`/`free`.
- **A fully container-native subtree:** `PPM_reconstruction_x/y`,
  `PPM_limit_pos`, `PPM_limit_cw84` and their `_fortran` halves.

**Scope — read this first.** In scope: changing a subroutine's array
dummies from `real, dimension(SZI_(G),…)` to `type(RealArray_t)` (or
`integer`/`logical` to `type(IntArray_t)`/`type(LogicalArray_t)`),
rewriting its body to obtain pointers via `%view`, and updating call
sites to build/write-back/free the containers. Out of scope: renaming
to `*_fortran`, a dispatcher shim, `select case (mode)`/`getenv_mode`,
capture mode/`io_recorder`, `#ifdef _TIM`, `bind(C)`, `%to_c`, anything
C++ — those belong to `generate_cpp_bridge`/`generate_amrex_code`,
layered on afterwards (§10). **A conversion never renames a
subroutine** — the `_fortran` examples above show what a container-
native subroutine looks like, not a naming convention this skill applies.

---

## 1. The big picture — moving the boundary, one routine at a time

Above the boundary, routines pass native Fortran arrays shaped by the
grid macros (`SZI_(G)`, `SZJ_(G)`, `SZK_(GV)`, …); below it, everything
speaks containers. Converting `FOO` moves the boundary up past it:

1. **`FOO` keeps its name.** Its array dummies become
   `type(RealArray_t)`/`type(IntArray_t)`/`type(LogicalArray_t)`, and
   its body obtains ordinary Fortran pointers via `%view`. Math untouched.
2. **Every caller takes on the marshalling** `FOO` used to get for
   free: `alloc` a container per array argument, call `FOO`, `copy2F`
   outputs back, `free` everything.

The marshalling doesn't disappear, it moves up a level; convert the
caller next and it moves again, until the boundary rests at the
module's genuine public entry points (in `MOM_continuity_PPM.F90`,
`continuity_PPM` and friends — no `Box_t`, probably meant to stay raw).

**Direction of travel.** One subroutine per invocation; nothing
recurses. **Top-down** (preferred) — start at the highest routine you
intend to convert; its caller is the one that stays raw permanently, so
the marshalling block is written once, and still-raw callees just take
a `%view` pointer in the meantime (§9 #10 shows the same trick for
`elemental` kernels). **Bottom-up** — leaves first; correct, but every
level's scaffolding gets written and then thrown away when that caller
is converted in turn. Either way, a converted routine calling a
still-raw callee is a valid, buildable intermediate state.

**This preference is about marshalling churn only — it doesn't decide
whether `G`/`GV`/`US` can be dropped from a signature.** Containerizing
a subroutine's own dummies is local; dropping a derived type
(`convert_array_containers` "Settle these decisions" #2, Step 2b) is
not — a subroutine that still forwards `G`/`GV` wholesale to even one
still-raw callee can't drop it, no matter how many of its own field
references have been promoted, and that's only knowable once every
callee below has finished the same check. So: containerize dummies
top-down, but decide `G`/`GV` elimination bottom-up, leaves first —
they're different questions, not a contradiction.

Because no signature *semantics* change, a conversion is **numerically
inert**; any numerical difference afterward is a bug.

---

## 2. Template — the converted subroutine (the callee)

Before:

```fortran
subroutine zonal_edge_thickness(bxC, h_in, h_W, h_E, G, GV, US, CS, OBC)
  type(box_t),             intent(in)    :: bxC
  type(ocean_grid_type),   intent(in)    :: G
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_in !< Tracer cell layer thickness [H ~> m or kg m-2].
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(out)   :: h_W  !< Western edge layer thickness [H ~> m or kg m-2].
  ...
  do concurrent (k=..., j=..., i=...)
    h_W(i,j,k) = h_in(i,j,k)
  enddo
```

After:

```fortran
subroutine zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                                h_min, upwind_1st, monotonic, simple_2nd, OBC)
  type(Box_t),          intent(in)    :: bxC        !< Iteration box for continuity solver
  type(RealArray_t),    intent(in)    :: h_in_a     !< Tracer cell layer thickness [H ~> m or kg m-2]
  type(RealArray_t),    intent(inout) :: h_W_a      !< Western edge layer thickness [H ~> m or kg m-2]
  type(RealArray_t),    intent(inout) :: h_E_a      !< Eastern edge layer thickness [H ~> m or kg m-2]
  type(RealArray_t),    intent(in)    :: mask2dT_a  !< Cell land/ocean mask [nondim]
  real,                 intent(in)    :: h_min      !< Minimum layer thickness (2*Angstrom_H) [H ~> m or kg m-2]
  logical,              intent(in)    :: upwind_1st !< If true, use 1st-order upwind reconstruction
  type(ocean_OBC_type), pointer       :: OBC        !< Open boundaries control structure

  integer :: i, j, k
  type(Box_t) :: bx
  real, dimension(:,:,:), contiguous, pointer :: h_in, h_W, h_E

  bx = bxC%grow(dim=1, n=1)
  call h_in_a%view(h_in)
  call h_W_a%view(h_W)
  call h_E_a%view(h_E)
  do concurrent (k=bx%idxS(3):bx%idxE(3), j=bx%idxS(2):bx%idxE(2), i=bx%idxS(1):bx%idxE(1))
    h_W(i,j,k) = h_in(i,j,k) ; h_E(i,j,k) = h_in(i,j,k)
  enddo
  call bx%free()
end subroutine zonal_edge_thickness
```

Key points:

- **Name unchanged.** No `_fortran`, no wrapper.
- **Loop body unchanged.** `%view` restores the container's `lb:ub`
  onto the pointer, so `h_W(i,j,k)` means exactly what it meant before
  — an edit inside the math means something is wrong.
- Container dummies are named `<orig>_a`; the `%view` pointer keeps the
  **original** name, which is what makes the body edit-free (§7).
- Scalars and `pointer` dummies (`OBC`) pass through unchanged. A
  grid-shaped `logical` *array* containerizes like `real`/`integer`, as
  `type(LogicalArray_t)`; only a bare `logical` *scalar* stays plain.
- `G`/`GV`/`US`/`CS` frequently disappear, replaced by the specific
  arrays/scalars the body actually used (§8).
- Written containers are `intent(inout)`, **never** `intent(out)` (§9 #6).

---

## 3. Template — the call site (the caller's new responsibility)

Every caller of a converted routine grows an alloc/call/write-back/free
block. Strip the dispatch machinery out of `zonal_edge_thickness`
(~636) and this is what remains:

```fortran
  real :: h_min
  type(RealArray_t) :: h_in_a, h_W_a, h_E_a, mask2dT_a

  h_min = 2.0 * GV%Angstrom_H

  ! Inputs: allocate AND copy in.
  call h_in_a%alloc(lb=LBOUND(h_in), ub=UBOUND(h_in), source=h_in)
  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  ! Pure outputs: copy in too -- the write below only covers the Box_t's
  ! range, not h_W/h_E's full halo-inclusive extent (§6).
  call h_W_a%alloc(lb=LBOUND(h_W), ub=UBOUND(h_W), source=h_W)
  call h_E_a%alloc(lb=LBOUND(h_E), ub=UBOUND(h_E), source=h_E)

  call zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                            h_min, CS%upwind_1st, CS%monotonic, CS%simple_2nd, OBC)

  ! Outputs and inouts: copy back out.
  call h_W_a%copy2F(h_W)
  call h_E_a%copy2F(h_E)

  call h_in_a%free()
  call h_W_a%free()
  call h_E_a%free()
  call mask2dT_a%free()
```

Fixed ordering: **derive scalars → alloc → call → `copy2F` → free.**
Leave any existing `cpu_clock_begin`/`cpu_clock_end` boundary where it is.

**This block is only for a caller that still takes raw arrays.** An
already container-based caller passes its own containers straight
through, no marshalling:

```fortran
  ! Caller already holds h_in_a / h_W_a / h_E_a as its own dummies.
  call zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                            h_min, upwind_1st, monotonic, simple_2nd, OBC)
```

Only an array the caller has no container for (a local work array,
say) still needs its own `alloc`/`free`.

So the block above is either **permanent** (the caller is a module
entry point staying raw) or **temporary scaffolding** (the caller is
itself queued for conversion, at which point the block is deleted).
Converting top-down avoids ever writing the temporary kind (§1).

> `h_W`/`h_E` are `intent(out)` yet the template still copies them in
> with `source=` — deliberate, not an oversight; see §6/§9 #2 for why
> "pure output" doesn't mean "safe to skip `source=`" here.

---

## 4. Container API — exact signatures

`RealArray_t`/`IntArray_t`/`LogicalArray_t` live in
`config_src/infra/FMS2/array_mod.F90` and
`config_src/infra/TIM/array_mod.F90` — the two infra layers are
**argument-for-argument identical** for everything here. `LogicalArray_t`
mirrors the other two with `logical` substituted for `real`/`integer`
— every call form, rule, and pitfall below (§4.1–§4.6) applies unchanged
(`allocLogical`, `copy2ALogical0D`–`4D`, `copy2FLogical1D`–`4D`, etc.).

### 4.1 `alloc` (generic)

Five specifics, disambiguated by the **rank of `source`**:

```fortran
! Scalar / no-source form -- `source` is OPTIONAL and SCALAR
subroutine allocReal(this, dims, lb, ub, source)
  class(RealArray_t), intent(inout) :: this
  integer, intent(in), optional :: dims(:)
  integer, intent(in), optional :: lb(:)
  integer, intent(in), optional :: ub(:)
  real(kind=real64), intent(in), optional :: source   ! fills every element

! Rank-N forms (1D/2D/3D/4D) -- `source` is REQUIRED and is an ARRAY
subroutine allocReal3D(this, dims, lb, ub, source)
  class(RealArray_t), intent(inout) :: this
  integer, intent(in), optional :: dims(:)
  integer, intent(in), optional :: lb(:)
  integer, intent(in), optional :: ub(:)
  real(kind=real64), intent(in) :: source(:,:,:)      ! deep copy-in
```

Four usable call forms:

| Call | Effect |
|---|---|
| `call a%alloc(lb=LBOUND(x), ub=UBOUND(x), source=x)` | allocate to `x`'s bounds **and deep-copy `x` in** |
| `call a%alloc(lb=LBOUND(x), ub=UBOUND(x))` | allocate only; contents **undefined** |
| `call a%alloc(lb=LBOUND(x), ub=UBOUND(x), source=0.0)` | allocate and fill with a scalar |
| `call a%alloc(dims=[n1,n2,n3])` | allocate with `lb=1`, `ub=dims` |

Runtime rules: supply **either** `lb`+`ub` **or** `dims`, never both,
never neither (else `FATAL "allocReal: Must specify either ub and lb or
dims"`); `size(lb) /= size(ub)` FATALs; `alloc` **deallocates any
existing payload first**, so re-allocating a live container is safe.

The array-`source` forms are what an older API spelled `%dup` followed
by `%copy2Array`; folded into `alloc` now (§4.6).

### 4.2 `allocView` (generic) — allocate and view in one call

```fortran
subroutine allocViewReal3D(this, a, dims, lb, ub, source)
   class(RealArray_t), intent(inout) :: this
   real(kind=real64), intent(inout), pointer :: a(:,:,:)   ! 2nd POSITIONAL arg
   integer, intent(in), optional :: dims(:)
   integer, intent(in), optional :: lb(:)
   integer, intent(in), optional :: ub(:)
   real(kind=real64), intent(in), optional :: source       ! SCALAR, optional
```

The pointer is the **second positional argument**, before
`dims`/`lb`/`ub`, and `source` here is scalar-and-optional even on the
rank-N specifics — the cleanest primitive for a pure output: allocate
and get a writable pointer without reading the caller's undefined array.

### 4.3 `view` (generic) — zero-copy pointer remap

```fortran
subroutine viewReal3D(this, a)
   class(RealArray_t), intent(in) :: this
   real(kind=real64), intent(inout), pointer :: a(:,:,:)
```

Body is literally `a(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)) => this%data`
— no allocation, no data movement. The pointer keeps the container's
original `lb:ub`, so ported index expressions are unchanged. `this` is
`intent(in)`, yet `view` still yields a **writable** pointer — the
in-tree code relies on this. Overloaded on the **rank of `a`** (1–4-D,
real and integer); a rank mismatch against `this%rank` FATALs at run
time. `a` is not declared `CONTIGUOUS`; passing a `contiguous, pointer`
actual is legal and standard here.

### 4.4 `copy2F` and `copy2Array` — directions

```fortran
subroutine copy2FReal3D(this, var)          ! container --> native Fortran array
  class(RealArray_t), intent(in) :: this
  real, dimension(:,:,:), intent(inout) :: var

subroutine copy2AReal3D(this, var)          ! native Fortran array --> container
  class(RealArray_t), intent(inout) :: this
  real, dimension(:,:,:), intent(in) :: var
```

Mnemonic: `copy2F` = copy *to Fortran* (out); `copy2Array` = copy *to
Array container* (in). `copy2F` exists for ranks 1–4 only (no 0-D);
`copy2Array` exists for ranks 0–4 (0-D broadcasts a scalar) and does
**not** allocate — the container must already exist, or it dereferences
a null `shape`.

### 4.5 `free`

```fortran
subroutine freeReal(this)
  class(RealArray_t), intent(inout) :: this
```

Deallocates payload plus all three metadata arrays, resets `rank = 0`.
Every deallocate is guarded by `associated`, so `free` is safe on an
unallocated or already-freed container, and calling it twice is safe.

### 4.5a `rebound`

```fortran
subroutine reboundReal(this, lb, ub)      ! also reboundInt, reboundLogical
  class(RealArray_t), intent(inout) :: this
  integer,            intent(in)    :: lb(:), ub(:)
```

Resets `lb`/`ub`/`shape` of an allocated container without reallocating or moving its data.
- It is FATAL if the container isn't allocated, if the rank differs, or if the new extents hold
  more elements than were allocated. Shrinking is allowed.
- Re-`%view` after a rebound: existing pointers keep their old bounds.
- **Use:** per-tile scratch whose shape is fixed but whose k-bounds move with the tile.
  Allocate once for the largest tile, then rebound and view inside the tile loop, as
  `CorAdCalc_TR` does. This replaces an alloc/free per tile.
- `size(this%data)` can exceed `product(shape)` after a shrink. `copy2F`/`copy2Array`, `%view`
  and the logical `to_c`/`from_c` all use `shape`, but `write_binary` writes all of `data`. So
  don't capture a shrunk container.

### 4.6 `dup` no longer exists

`%dup` was restructured into `%alloc`: `call a%dup(x) ; call
a%copy2Array(x)` is now the single `call a%alloc(lb=LBOUND(x),
ub=UBOUND(x), source=x)`. No `dup` type-bound procedure exists in
either infra layer (`grep -rn "%dup"` returns zero hits) — anything
still showing `%dup` predates this and won't compile.

---

## 5. `Box_t` — iteration domain

Defined in `config_src/infra/{FMS2,TIM}/box_mod.F90`. The two files are the
same except for where `Box_C` comes from: FMS2's defines it, and TIM's
imports it from TIM's `turbotmp_bridge_c_types` module and re-exports it.
The `Box_t` API below is identical in both.
`idxS(:)`/`idxE(:)` are `allocatable` and public; every binding is a
plain `procedure ::` (no generics, one signature each).

| Call | Semantics |
|---|---|
| `call bx%safe_alloc(ndims)` | Deallocate if needed, then allocate `idxS`/`idxE` to `ndims`, zero-filled. Idempotent. |
| `call bx%set(idxS=[...], idxE=[...])` | Element-wise assign. **FATALs if not already allocated** — `safe_alloc` must come first. |
| `call bx%free()` | Deallocate both components. Safe on an unallocated box. |
| `new = bx%grow(dim, n)` | **Returns a NEW `Box_t` by value.** `idxS(dim) -= n`, `idxE(dim) += n`. Does not mutate `bx`. |
| `new = bx%growLo(dim, n)` / `growHi(dim, n)` | As `grow`, but only the start / only the end. |
| `new = bx%growBy(lo=[…], hi=[…])` | **Returns a NEW `Box_t`.** `idxS(:) -= lo(:)`, `idxE(:) += hi(:)`, one entry per dimension; FATAL if the sizes don't match the rank. Use it for any multi-dimension or asymmetric grow, instead of chaining `grow*` through a temporary. In AMReX this is `Box(bx.smallEnd() - lo, bx.bigEnd() + hi)`, with dimensions numbered from 0. |
| `new = bx%shrink(dim, n)` | New box with `idxS(dim) += n`, `idxE(dim) -= n`. |
| `cdesc = bx%to_c()` | `Box_C` of `c_ptr`s. Null-safe. Available under **both** infra layers. |

`dim` is 1-based, not bounds-checked. Axis convention: **`dim=1` is
zonal/i/x, `dim=2` is meridional/j/y** — the key difference between a
zonal routine and its meridional mirror.

Standard iteration idiom:

```fortran
do concurrent (k=bx%idxS(3):bx%idxE(3), j=bx%idxS(2):bx%idxE(2), i=bx%idxS(1):bx%idxE(1))
```

`grow`/`shrink` return a **fresh allocation you own**:

```fortran
bx = bxC%grow(dim=1, n=1)
...
call bx%free()          ! required -- see §9 #3
```

**Reading bounds once at the top of a kernel.** Instead of growing a
new box for each loop, a kernel can read scalar bounds from its boxes
once, the way the continuity kernels do
(`ish = bxC%idxS(1) ; ieh = bxC%idxE(1) ; …`), and keep its loop
headers' original range expressions (`I=ish-1:ieh`). That is the
lower-risk choice when converting an existing body: every loop header
can be checked against the original line by line. Once the tree is
stable, `convert_loops_to_box_iterators` replaces those scalar ranges
with core and derived boxes (e.g. `bxQ_pij = bxQ` grown at the end in i
and j) ahead of AMReX bridging.

```fortran
is  = bxH%idxS(1) ; ie  = bxH%idxE(1) ; js  = bxH%idxS(2) ; je  = bxH%idxE(2)
Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
ksc = bxH%idxS(3) ; kec = bxH%idxE(3)      ! not ks/ke: see §9 #18
```

**Memory mode: build staggered boxes from `G`, not by growing.**
`G%IscB`/`G%JscB` equal `isc`/`jsc` with non-symmetric memory and
`isc-1`/`jsc-1` with symmetric memory. So growing an h-point box by one
reproduces a B-index range (`Isq:Ieq`) in only one of the two modes.
When the code being converted uses `G%IscB`-style ranges, build a
separate box from them, e.g.
`bxQ%set(idxS=[G%IscB,G%JscB,ks], idxE=[G%IecB,G%JecB,ke])`, so both
memory modes stay bit-for-bit. (Continuity avoids the issue because its
code always used `ish-1`-style offsets.)

**Per-tile boxes and tile-shaped scratch.** A root routine that loops
over tiles (k-blocks today, AMReX `MFIter` tiles later) allocates its
boxes once and re-`set`s them per tile:

```fortran
call bxH%safe_alloc(ndims=3)
do k_start = 1, nz, nkblock
  k_end = min(k_start+nkblock-1, nz)
  call bxH%set(idxS=[G%isc,G%jsc,k_start], idxE=[G%iec,G%jec,k_end])
  call tile_body(..., bxH, ...)
enddo
call bxH%free()
```

Inside the per-tile body, scratch can be an automatic array whose k-bounds
come from the box, indexed by the real `k` rather than a block-relative
`kk`. This is legal because a specification expression may reference a
component of a non-optional, non-`intent(out)` dummy:

```fortran
type(Box_t), intent(in) :: bxH    ! declare the box before the arrays that use it
real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)) :: q
```

The per-tile body's writes to full-k outputs (`CAu`, `CAv`) cover only
that tile's layers, so those dummies must be `intent(inout)`, not
`intent(out)`: an `intent(out)` array would lose the other tiles'
layers.

**Persistent grid containers: `grid_core` and `grid_OBC`** (`src/core/MOM_grid_containers.F90`).

- **What they are:**
  - Copies of the `ocean_grid_type` arrays, built once in `initialize_MOM` (after the grid is
    final) and freed in `MOM_end`.
  - `grid_core` holds every non-OBC array a planned tree reads.
  - `grid_OBC` holds `OBCmaskCu IdxCu_OBCmask OBCmaskCv IdyCv_OBCmask`.
  - Members are `RealArray_t` with the **same names as in `G`**. Each structure also carries `C`,
    its `bind(C)` mirror (`grid_core_C` / `grid_OBC_C`), built once under `#ifdef _TIM`.
- **Getting them:** `Gcore => grid_core(G)` / `Gobc => grid_OBC(G)` return pointers. They are
  FATAL if the structure was not built, or if its bounds don't match `G`.
- **Using them in a kernel:**

  ```fortran
  type(grid_core_type), intent(in) :: Gcore !< Persistent copies of the core ocean grid fields
  real, dimension(:,:), contiguous, pointer :: IdxCu, IdyCv
  call Gcore%IdxCu%view(IdxCu) ; call Gcore%IdyCv%view(IdyCv)
  ```

  The body then uses `IdxCu` where it had `G%IdxCu`.
- **Never** `%free()` a member, and never `alloc(source=G%<f>)` per call.
- **Adding members:** new members are appended, in the same order in all five places (both
  types, init, `to_c`, end). Their comments are copied word for word from `MOM_grid.F90`. Adding
  a field, or a new group, is the user's decision.
- **Tree-root boxes:** a tree root gets grid index ranges as whole-column boxes from its cap
  (`bxH0`, `bxQ0`), and builds U/V tile boxes by combining their i/j ranges.

**On the C++ side** every array's own `lb` positions its `Array4`. See
protoMOMxx `generate_amrex_code` lessons.md §7 #12: arrays with different
staggering in one kernel break silently if that's skipped.

---

## 6. Intent-driven copy rules

What a dummy's original `intent` implies for the call-site block (§3):

| Original intent | Allocate | Copy in | Copy back (`copy2F`) | Container dummy intent |
|---|---|---|---|---|
| `intent(in)` | yes | **yes** — `source=x` | no | `intent(in)` |
| `intent(out)` | yes | **default yes** — `source=x` (see below) | **yes** | `intent(inout)` |
| `intent(inout)` | yes | **yes** — `source=x` | **yes** | `intent(inout)` |

**`intent(out)` does not mean "the callee writes every element."** It
only means the array has no meaningful value on entry, and for a plain
array it does **not** clear memory — before conversion, whatever the
callee didn't write, the caller kept from before.

**This matters here because almost every array is allocated wider than
any single call writes it.** `SZI_`/`SZIB_`/`SZJ_`/`SZK_`-shaped arrays
span the full halo-inclusive data domain, while a `Box_t`-scoped write
covers only the narrower computational domain — the gap is the halo,
untouched by most calls on purpose. Allocate without `source=` for such
a write and the halo starts as garbage; `copy2F` then copies that
garbage back over the caller's real array, silently destroying whatever
meaningful halo content was there. Not hypothetical: this is the root
cause of a real, hard-to-trace CI failure in an already-merged routine
(`zonal_mass_flux`'s `uh`) — it doesn't crash, the corrupted halo cells
just sit wrong until something downstream reads them.

**Default: copy in with `source=x`, even for `intent(out)`, unless you
can positively confirm the callee writes the array's full `LBOUND`/
`UBOUND`** (not just the `Box_t`'s range) — for a `Box_t`-scoped write
that confirmation essentially never holds here. The allocate-only forms
below are real but are the exception:

```fortran
call h_W_a%alloc(lb=LBOUND(h_W), ub=UBOUND(h_W))                  ! only if h_W is written in full
call h_W_a%allocView(h_W_p, lb=LBOUND(h_W), ub=UBOUND(h_W))       ! same caveat -- no copy-in either
```

If in doubt, copy in — a wasted halo copy is cheap, a silently
corrupted one is not.

### Optional dummies

An `optional` array dummy keeps `optional` when it becomes a container,
and the `present()` machinery moves with it.

**In the callee**, three things change together:

```fortran
  type(RealArray_t), optional, intent(in) :: hin_a  !< Initial layer thickness [H ~> m or kg m-2].
                                                    !! If hin is absent, h is also the initial thickness.
  real, dimension(:,:,:), contiguous, pointer :: hin

  if (present(hin_a)) call hin_a%view(hin)      ! guard the view
  ...
  if (present(hin_a)) then                      ! existing logic now tests the _a name
```

1. Keep `optional` on the container declaration.
2. **Guard the `%view`** — calling it on an absent dummy is illegal.
3. **Rewrite every `present(<orig>)` test to `present(<orig>_a)`** —
   easy to miss, since a stale `present(hin)` won't compile once `hin`
   is a local pointer rather than a dummy.

The pointer stays unassociated when absent, so any use must sit inside
the same `present` guard.

**At the call site**, it depends on whether the caller's own array is
optional:

- **Always present** (the common case): allocate and pass normally.
- **Caller's array is itself `optional`**: a container local is
  *always* present as an actual argument, so passing an unallocated one
  would make `present()` true in the callee with garbage behind it.
  Allocate conditionally and branch the call:

  ```fortran
  if (present(hin)) then
    call hin_a%alloc(lb=LBOUND(hin), ub=UBOUND(hin), source=hin)
    call continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a, hin_a=hin_a)
  else
    call continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a)
  endif
  call hin_a%free()      ! safe even if never allocated (§4.5)
  ```

  Fine with one or two optionals; past two, the branch count explodes —
  **stop and ask the user** rather than nesting conditionals.

Optional *scalars* (`hmin`, `h_min`) are unaffected — never
containerized, `optional`/`present()` unchanged.

---

## 7. Naming conventions

| Thing | Convention | Example |
|---|---|---|
| Subroutine name | **unchanged** | `zonal_mass_flux` stays `zonal_mass_flux` |
| Container dummy | original name + `_a` | `h_in` → `h_in_a` |
| `%view` pointer inside the body | the **original** name | `h_in` |
| Container local at a call site | original name + `_a` | `h_W_a` |
| Grid-derived container | field name + `_a` | `G%mask2dT` → `mask2dT_a` |

Giving the `%view` pointer the original name is what keeps the loop
body edit-free — a clean conversion diff touches declarations and adds
`%view` calls, nothing in the math.

### Doc comments

Every dummy carries a `!<` Doxygen comment, and **every one must
survive a conversion verbatim** — same wording, unit annotation, and
trailing punctuation. The common way this breaks: the original spans
two lines with the comment on the continuation line, while the
converted declaration fits on one, and the comment gets lost in the
collapse:

```fortran
! BEFORE -- two lines, comment on the continuation
  real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                           intent(in)    :: h_in !< Tracer cell layer thickness [H ~> m or kg m-2].

! AFTER -- one line, comment carried across
  type(RealArray_t),    intent(in)    :: h_in_a     !< Tracer cell layer thickness [H ~> m or kg m-2]
```

Record each comment before rewriting the declaration, then re-attach
it. A dummy **added** by the conversion (a grid-derived container, a
scalar lifted out of `CS`/`GV`) needs a new `!<` comment in house
style — a missing one is a Doxygen build warning, and CI treats doc
coverage as part of the build. Comments not attached to a dummy (the
`!>` subroutine header, blocks between declarations, body comments)
stay exactly as they are.

### No commentary about the conversion itself

This codebase documents ocean physics, never the fact that a subroutine
was migrated or why a Fortran mechanism is used here. Don't add
comments like `! Containers for zonal_mass_flux` or `! Must be
nullified -- see below` anywhere — the API and `_a` convention already
make the pattern legible; a conversion diff should read as pure
mechanism, not narration. Exception, already covered above: a
genuinely new dummy still needs its own `!<` doc comment describing its
physical meaning.

---

## 8. Grid-derived arrays and shrinking argument lists

> **`ocean_grid_type` arrays are the exception to this section.** They come from the persistent
> `grid_core`/`grid_OBC` containers (§5, "Persistent grid containers"), not from per-field
> dummies built by the caller. What follows still applies to other derived types (`GV`, `US`,
> `CS`, …).

Some arrays a routine needs aren't dummies at all — reached through a
derived type, most commonly `G%mask2dT`. These become **new** container
dummies, built by the caller:

```fortran
call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
```

So a converted argument list isn't the old list with types swapped —
it gains container dummies and often sheds `G`/`GV`/`US`/`CS`. Classify
every derived-type reference in the original body:

- **Array** (`G%mask2dT`, `G%areaT`) → new container dummy, allocated
  by the caller.
- **Scalar** (`GV%Angstrom_H`, `CS%upwind_1st`) → plain scalar dummy,
  computed at the call site (`h_min = 2.0 * GV%Angstrom_H`).
- **Still genuinely needed whole** (forwarded to an unconverted callee,
  or resists splitting) → keep the derived type. Don't force it.

Only drop `G`/`GV`/`US`/`CS` once nothing in the body references them.

---

## 9. Recurring pitfalls

1. **Do not rename the subroutine.** No `_fortran` suffix, no wrapper
   pair — that's `generate_cpp_bridge`'s job (§10).
2. **Skipping `source=` on a `Box_t`-scoped "pure output" silently
   corrupts its halo.** `intent(out)` only means undefined on entry, not
   written in full; a `Box_t`-scoped write covers a narrower range than
   the array's full halo-inclusive extent, so allocate-without-`source=`
   leaves the halo as garbage that `copy2F` then copies back over the
   caller's real array. Caused a real, hard-to-trace CI failure
   (`zonal_mass_flux`'s `uh`) — default to `source=x` even for
   `intent(out)` (§6); skip the copy-in only once confirmed the write
   covers the full `LBOUND`/`UBOUND`, not just the `Box_t`'s range.
3. **`grow`/`shrink`/`growLo`/`growHi` return a NEW box you must
   `free`.** Forgetting leaks on every call — inside timestep loops.
4. **`set` before `safe_alloc` is FATAL.** Always `safe_alloc(ndims=3)`
   first.
5. **`copy2Array` does not allocate.** `alloc(..., source=x)` is the
   one-step form; `copy2Array` only re-fills an existing container.
6. **Written container dummies must be `intent(inout)`, never
   `intent(out)`.** `intent(out)` resets the pointer components,
   discarding the caller's allocation — inverts the raw-array
   convention and is easy to get backwards.
7. **Container extents must match the source array's extents.** Copy
   loops are driven by `this%shape`, not `SHAPE(source)`; a mismatch
   reads/writes out of bounds instead of failing cleanly.
8. **`view` FATALs on rank mismatch at run time, not compile time.**
8a. **`optional` dummies: rename the `present()` tests too.** A missed
    rename doesn't compile once `hin` is a local pointer; a missed
    *guard* compiles and fails at run time.
8b. **Never pass an unallocated container to an `optional` dummy** — a
    container local is always present as an actual argument, so
    `present()` reads true regardless. Branch the call if the caller's
    source is itself optional (§6).
9. **Every call site must be updated in the same change** — find them
   with `grep -rn "call[[:space:]]\+<name>(" src/ config_src/`,
   checking continuation lines too.
10. **`elemental` callees cannot take containers.** `flux_elem`/
    `flux_elem_OBC` keep their scalar signatures; a converted caller
    passes `%view` pointer elements to them, unchanged.
11. **Derived types with array members are not containers.**
    `BT_cont_type` (`%h_u`, `%h_v`) holds raw arrays; converting a
    routine that passes `BT_cont` doesn't convert its members — leave
    the type alone unless explicitly asked.
12. **`%to_c` on the array types is TIM-only** (FMS2 has it commented
    out) — every array `%to_c()` must sit inside `#ifdef _TIM`
    (`Box_t%to_c` is available under both). A conversion emits no
    `to_c` at all; this only matters once a bridge is added later.
13. **`IntArray_C` is not public** in either infra layer, so a caller
    can't declare the result type of `to_c_Int`.
14. **Default `real` is assumed `real64`.** `copy2*` declares `var` as
    default `real` while `alloc*` declares `source` as
    `real(kind=real64)` — compiles only because MOM6 builds with
    default real = real64. Never mix kinds deliberately across this.
15. **Free every container on every path** — early returns must still
    reach the `free` block.
16. **Every pointer local that is only conditionally `%view`'d must be
    `nullify`'d before the first guard that might set it** — whether or
    not it's later forwarded, whether or not every current read sits in
    the same guard. A plain local pointer's association status on entry
    is **undefined, not disassociated**; skip the `nullify` and it's
    left holding garbage, with no compile error and usually no crash at
    the point of the mistake — confirmed as a real, CI-caught silent
    numerical-corruption bug here
    (`present_uhbt_or_set_BT_cont`/`present_vhbt_or_set_BT_cont`). Do
    **not** nullify in the declaration (`=> null()`) — that implicitly
    adds `SAVE`, so it only runs once ever.

    One sub-variety: forwarding the pointer *onward* to a still-raw
    callee's `optional` dummy. A `optional, type(RealArray_t)` can't be
    passed directly to a raw optional dummy (type mismatch); fix is the
    same guarded `%view`, forwarded *unconditionally* by keyword — a
    disassociated pointer reads as genuinely absent to a non-pointer,
    non-allocatable `optional` dummy, but only if nullified first, and
    only if the callee's own dummy is a plain optional array (never
    pointer/allocatable). See `convert_optional_args_to_containers` for
    the full mechanism.
17. **Never reference a container's value — `%associated()` or any
    other method/component — from *inside* an `!$omp
    target`/`teams`/`parallel`/`loop`-family construct.** Precompute a
    plain scalar before the construct and reference only that inside.
    Distinct from #16 (uninitialized *pointer* vs. the *container*
    itself). Root-caused via AddressSanitizer:
    `zonal_flux_adjust_fortran`/`meridional_flux_adjust_fortran` called
    `uh_3d_a%associated()`/`vh_3d_a%associated()` directly inside a
    `!$omp target teams loop`; macOS has no real GPU so `target` falls
    back to `libgomp`, but a variable referenced inside the construct
    without an explicit `map`/`private`/`firstprivate` still gets
    implicitly captured, and a container's pointer components confused
    the outlined region's cleanup — ASan reported a `bad-free` landing
    inside an unrelated caller-local's stack slot. Not specific to a
    never-allocated container; the trigger is referencing the container
    inside the construct at all. Fix: the pattern this same file
    already used correctly for `uhbt_a` — `use_uhbt =
    uhbt_a%associated()` computed once outside the loop, only
    `use_uhbt` referenced inside. Audit for this the same way as a
    missing `nullify` (#16) — both are "prepare the value entirely
    before the construct begins."
18. **Fortran names are case-insensitive: check every new name against
    existing symbols ignoring case.** A new `ke` (a tile's k-end bound)
    is the same identifier as an existing array `KE` (kinetic energy).
    In CorAdCalc's `CorAdv_tile` and `gradKE` this failed to compile,
    and the error was confusing:
    - the real error is `Symbol 'ke' at (1) already has basic type of
      REAL`, at the declaration;
    - gfortran also reports cascade errors at unrelated lines, e.g.
      `Cannot convert INTEGER(4) to UNKNOWN` at a
      `do concurrent (k=ks:ke, …)` header.

    Before introducing a dummy, local or loop bound, grep the subroutine
    for it with `grep -niw`, and run a case-insensitive duplicate-
    declaration check per edited subroutine before handing the code
    over. When a build fails with several errors, look for the first
    one: later errors are often cascades. Names that avoid this here:
    `ksc`/`kec` (like MOM6's `isc`/`iec`).

---

## 10. Verification, and composing with the other skills

A conversion is **numerically inert by construction** (§1), so
verification is mostly about proving nothing moved:

1. **Diff review.** Math unchanged; any edit inside a loop body is a
   red flag. The diff is dummy declarations, added `%view` calls, and
   call-site alloc/copy2F/free blocks.
2. **Call-site completeness.** Every caller found by grep is updated;
   argument count and order match the new signature.
3. **Both infra layers.** Must compile under FMS2 *and* TIM — the API
   is identical across them (§4), so a single-layer-only compile means
   something reached outside the shared surface.
4. **Build — the user's job.** Don't compile, build or run tests
   yourself, not even a single-file syntax check, unless the user
   explicitly allows it in this conversation (user rule,
   `~/.claude/CLAUDE.md`). Say plainly that the build is pending and
   what the user should build and run. Never call a conversion
   verified on inspection alone.
5. **Numerics.** Bit-identical to the pre-conversion build.

**Composing with the bridge skills.** A converted routine is the
natural input to `generate_cpp_bridge`, which adds the `*_fortran`
rename, dispatcher shim, `getenv_mode`, capture mode, the `#ifdef _TIM`
AMReX arm, and the `bind(C)` interface; `generate_amrex_code` then
writes the C++ side. Keeping the conversion separate keeps each diff
small and reviewable, and a routine can be containerized with no
commitment to ever bridging it.

**A `_fortran` sibling doesn't change what this skill does.** A
subroutine with raw dummies is a conversion target regardless of
history — never touched, or a bridge shim whose signature hasn't been
migrated yet. Convert it exactly like any other raw-array subroutine
(§2–§8), with one constraint: some of its content belongs to
`generate_cpp_bridge`, not this skill, and must be left alone —

- Never edit anything between `select case (mode)` and `end select`,
  except renaming a `CS%x`/`GV%x` reference if that derived type is
  dropped (§8). Nothing else in that region changes.
- Never touch a `bind(C) interface … end interface` declaration — it
  takes `RealArray_C`/`Box_C` regardless of where the `RealArray_t`
  came from, so a conversion never needs to reach it.

If you find yourself rewriting either, stop — the task has drifted into
bridge territory. Container-building/writeback bookkeeping *around*
that content is ordinary conversion scaffolding and gets replaced the
same way as any call site (§3) — it just happens to be inside this
subroutine rather than a caller.

`meridional_edge_thickness` (not yet migrated) and `zonal_edge_thickness`
(migrated) are an exact real-world before/after pair, if you want to
see it worked out in full.

The reliable signal that a raw-array subroutine has bridge content to
preserve is **not** its dummy types — it's a `<name>_fortran` sibling's
existence and dispatch machinery in the body (`select case (mode)`,
`getenv_mode`, `TIMH_*`, `io_recorder`, `#ifdef _TIM`), which
`convert_array_containers`'s Step 0 already checks for.
