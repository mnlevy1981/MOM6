---
name: array_container_lessons
version: "0.4.0"
description: Reference for converting MOM6 subroutines from raw grid-shaped Fortran arrays to RealArray_t / IntArray_t / LogicalArray_t / Box_t containers -- callee and call-site templates, the exact alloc/allocView/view/copy2F/copy2Array/free/rebound and Box_t API, intent-driven copy rules, naming and doc-comment rules, and numbered pitfalls, in sections §1-§10 that convert_array_containers and its sibling skills cite by number. Invoke once per session before the first convert_array_containers run. STANDING RULE for every MOM6 edit: never reword, shorten or "fix" an existing Fortran comment; copy doc comments verbatim (see the box at the top).
user-invocable: true
---

# Array-container conversion: reference

> **STANDING RULE — Fortran comments are not yours to edit (user, repeated; 2026-09-29).**
> - Never change the wording of an existing Fortran comment (`!<`, `!!` or `!`): no shortening,
>   paraphrasing, dropped article or trailing period, or "corrected" typo.
> - A declaration whose doc comment comes from elsewhere (a `grid_core` member from
>   `MOM_grid.F90`, a promoted dummy from the field it replaces) copies it **word for word**, even
>   where it looks wrong (`IdyCv_OBCmask` says "1/dxCv").
> - Too long? Move line breaks only, continuing on `!!` lines aligned under the `!<`. Keep a
>   bracketed unit, and alternatives like `[degrees_N] or [km] or [m]`, on one line.
> - Before handing over, check by script that each copied comment equals its source word for word.
> - If a change makes an existing comment inaccurate, leave it and tell the user.
> - A new comment on a new declaration is fine; editing an existing one is not.

> **STANDING RULE — use the pre-written skills (user, 2026-10-23).** Every MOM6 conversion step a
> skill in `.claude/skills/` covers is done by invoking that skill (Skill tool) and following its
> procedure. A script may apply a skill's steps mechanically but never does anything the skill
> doesn't say. If a skill lacks something, propose a general addition for the user to review
> instead of improvising one-off code. A campaign is driven by invoking `convert_calltree`, not by
> reading its file; hand-authored code is allowed only where `convert_calltree` Step 5a says so
> (Stage 1).

Sections §1–§10 are cited by number from other skills; do not renumber them. Every touched file
carries exactly one `!!SKILLS: 0.4.0` marker after the license header.

Worked reference: `src/core/MOM_continuity_PPM.F90` — `zonal_edge_thickness_fortran` /
`meridional_edge_thickness_fortran` (converted callees), `zonal_edge_thickness` /
`meridional_edge_thickness` (call-site marshalling), `PPM_reconstruction_x/y`, `PPM_limit_pos`,
`PPM_limit_cw84` (a container-native subtree). The `_fortran` names come from a later bridge; a
conversion itself never renames.

**Scope:** array dummies become containers, the body takes `%view` pointers, call sites marshal.
Renaming, dispatcher shims, capture mode, `#ifdef _TIM`, `bind(C)`, `%to_c` and C++ belong to
`generate_cpp_bridge` (§10).

---

## 1. The big picture

Above the boundary, routines pass raw arrays sized by grid macros (`SZI_(G)`, `SZK_(GV)`, …);
below it, containers. Converting `FOO` moves the boundary up past it: `FOO` keeps its name, its
array dummies become containers, its body views them; every caller takes on the marshalling
(alloc → call → `copy2F` → free). Converting the caller moves it up again, until it rests at the
module's genuine public entry points.

- **Containerizing dummies may go top-down** (preferred: the permanently-raw caller's
  marshalling is written once; still-raw callees get a `%view` pointer meanwhile). Bottom-up is
  correct but writes scaffolding that is later thrown away. A converted routine calling a
  still-raw callee is a valid intermediate state.
- **Dropping `G`/`GV`/`US`/`CS` is decided bottom-up**: a routine can't drop a derived type while
  it still forwards it to any callee, and that is known only once every callee has settled.
- A conversion changes no semantics, so it is **numerically inert**; any difference is a bug.

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
  type(RealArray_t),    intent(in)    :: h_in_a     !< Tracer cell layer thickness [H ~> m or kg m-2].
  type(RealArray_t),    intent(inout) :: h_W_a      !< Western edge layer thickness [H ~> m or kg m-2].
  type(RealArray_t),    intent(inout) :: h_E_a      !< Eastern edge layer thickness [H ~> m or kg m-2].
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

- Name unchanged. Loop body unchanged: `%view` keeps the container's `lb:ub`, and the pointer
  takes the original name (§7). An edit inside the math means something is wrong.
- Grid-shaped `real`/`integer`/`logical` arrays → `RealArray_t`/`IntArray_t`/`LogicalArray_t`.
  Scalars (including `logical` scalars) and `pointer` dummies (`OBC`) pass through unchanged.
- Derived-type dummies are replaced by what the body actually used (§8).
- Written containers are `intent(inout)`, never `intent(out)` (§9 #6).

---

## 3. Template — the call site

A caller that still takes raw arrays (Case A) marshals:

```fortran
  real :: h_min
  type(RealArray_t) :: h_in_a, h_W_a, h_E_a, mask2dT_a

  h_min = 2.0 * GV%Angstrom_H

  call h_in_a%alloc(lb=LBOUND(h_in), ub=UBOUND(h_in), source=h_in)
  call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
  call h_W_a%alloc(lb=LBOUND(h_W), ub=UBOUND(h_W), source=h_W)   ! outputs copy in too (§6)
  call h_E_a%alloc(lb=LBOUND(h_E), ub=UBOUND(h_E), source=h_E)

  call zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, &
                            h_min, CS%upwind_1st, CS%monotonic, CS%simple_2nd, OBC)

  call h_W_a%copy2F(h_W)
  call h_E_a%copy2F(h_E)

  call h_in_a%free() ; call h_W_a%free() ; call h_E_a%free() ; call mask2dT_a%free()
```

Order: **derive scalars → alloc → call → `copy2F` → free.** Existing
`cpu_clock_begin`/`cpu_clock_end` calls stay where they are. (For `G` arrays a converted tree
views `grid_core` instead of allocating per call — §5, §8.)

A caller that is already converted (Case B) passes its own containers straight through; only an
array it holds no container for (a local work array) needs its own alloc/free. The Case A block is
therefore permanent (caller stays raw) or temporary scaffolding (caller is queued for conversion).

---

## 4. Container API — exact signatures

`RealArray_t`/`IntArray_t`/`LogicalArray_t` live in `config_src/infra/{FMS2,TIM}/array_mod.F90`;
the two layers are argument-for-argument identical for everything here. `LogicalArray_t` mirrors
the others with `logical` substituted (`allocLogical`, `copy2ALogical0D`–`4D`,
`copy2FLogical1D`–`4D`, …); every rule below applies to all three.

### 4.1 `alloc` (generic)

Five specifics, chosen by the **rank of `source`**:

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

| Call | Effect |
|---|---|
| `call a%alloc(lb=LBOUND(x), ub=UBOUND(x), source=x)` | allocate to `x`'s bounds **and deep-copy `x` in** |
| `call a%alloc(lb=LBOUND(x), ub=UBOUND(x))` | allocate only; contents **undefined** |
| `call a%alloc(lb=LBOUND(x), ub=UBOUND(x), source=0.0)` | allocate and fill with a scalar |
| `call a%alloc(dims=[n1,n2,n3])` | allocate with `lb=1`, `ub=dims` |

Supply `lb`+`ub` **or** `dims`, never both or neither (else `FATAL "allocReal: Must specify
either ub and lb or dims"`); `size(lb) /= size(ub)` FATALs. `alloc` deallocates any existing
payload first, so re-allocating a live container is safe.

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

The pointer is the second positional argument; `source` is scalar and optional even on the
rank-N specifics. Gives a writable pointer without reading the caller's array (but see §6).

### 4.3 `view` (generic) — zero-copy pointer remap

```fortran
subroutine viewReal3D(this, a)
   class(RealArray_t), intent(in) :: this
   real(kind=real64), intent(inout), pointer :: a(:,:,:)
```

Body is `a(lb(1):ub(1), lb(2):ub(2), lb(3):ub(3)) => this%data`: no allocation or copy, and the
pointer keeps the container's `lb:ub`. `this` is `intent(in)` yet the pointer is writable (in-tree
code relies on this). Overloaded on the rank of `a` (1–4-D); a mismatch with `this%rank` FATALs
at run time. `a` isn't `CONTIGUOUS`; passing a `contiguous, pointer` actual is legal.

### 4.4 `copy2F` and `copy2Array` — directions

```fortran
subroutine copy2FReal3D(this, var)          ! container --> native Fortran array
  class(RealArray_t), intent(in) :: this
  real, dimension(:,:,:), intent(inout) :: var

subroutine copy2AReal3D(this, var)          ! native Fortran array --> container
  class(RealArray_t), intent(inout) :: this
  real, dimension(:,:,:), intent(in) :: var
```

`copy2F` = to Fortran (out), ranks 1–4. `copy2Array` = to Array (in), ranks 0–4 (0-D broadcasts
a scalar); it does **not** allocate — on an unallocated container it dereferences a null `shape`.

### 4.5 `free`

```fortran
subroutine freeReal(this)
  class(RealArray_t), intent(inout) :: this
```

Deallocates the payload and metadata and resets `rank = 0`. Guarded by `associated`, so safe on an
unallocated or already-freed container.

### 4.5a `rebound`

```fortran
subroutine reboundReal(this, lb, ub)      ! also reboundInt, reboundLogical
  class(RealArray_t), intent(inout) :: this
  integer,            intent(in)    :: lb(:), ub(:)
```

Resets `lb`/`ub`/`shape` of an allocated container without moving data. FATAL if unallocated, if
the rank differs, or if the new extents hold more elements than allocated (shrinking is fine).
Re-`%view` afterwards; existing pointers keep the old bounds. Use: per-tile scratch allocated once
for the largest tile and rebounded per tile (as `CorAdCalc_TR` does). After a shrink
`size(data) > product(shape)`; `copy2*`, `%view`, logical `to_c`/`from_c` use `shape` but
`write_binary` writes all of `data`, so don't capture a shrunk container.

### 4.6 `dup` no longer exists

`call a%dup(x) ; call a%copy2Array(x)` is now `call a%alloc(lb=LBOUND(x), ub=UBOUND(x),
source=x)`. No `dup` binding exists in either layer; code using `%dup` won't compile.

---

## 5. `Box_t`, tiles and persistent grid containers

`Box_t` is in `config_src/infra/{FMS2,TIM}/box_mod.F90` (same API; TIM's imports `Box_C` from
TIM's `turbotmp_bridge_c_types`). `idxS(:)`/`idxE(:)` are public allocatables; each binding has
one signature.

| Call | Semantics |
|---|---|
| `call bx%safe_alloc(ndims)` | (Re)allocate `idxS`/`idxE` to `ndims`, zero-filled. Idempotent. |
| `call bx%set(idxS=[...], idxE=[...])` | Element-wise assign; **FATAL unless `safe_alloc` came first.** |
| `call bx%free()` | Deallocate; safe on an unallocated box. |
| `new = bx%grow(dim, n)` | **New box** with `idxS(dim) -= n`, `idxE(dim) += n`; `bx` unchanged. |
| `new = bx%growLo(dim, n)` / `growHi(dim, n)` | As `grow`, start only / end only. |
| `new = bx%growBy(lo=[…], hi=[…])` | **New box**, `idxS(:) -= lo(:)`, `idxE(:) += hi(:)`; FATAL if sizes ≠ rank. Use for any multi-dimension or asymmetric grow instead of chaining. AMReX: `Box(bx.smallEnd() - lo, bx.bigEnd() + hi)`, dims from 0. |
| `new = bx%shrink(dim, n)` | New box with `idxS(dim) += n`, `idxE(dim) -= n`. |
| `cdesc = bx%to_c()` | `Box_C` of `c_ptr`s; null-safe; both infra layers. |

`dim` is 1-based and not bounds-checked; `dim=1` is i/x (zonal), `dim=2` j/y (meridional). Every
returned box must be freed (§9 #3). Iteration idiom:
`do concurrent (k=bx%idxS(3):bx%idxE(3), j=bx%idxS(2):bx%idxE(2), i=bx%idxS(1):bx%idxE(1))`.

**Scalar bounds read once.** When converting an existing body, read scalar bounds from the boxes
at the top and keep the original loop-header expressions, so each header checks line by line;
`convert_loops_to_box_iterators` later turns them into box iterations.

```fortran
is  = bxH%idxS(1) ; ie  = bxH%idxE(1) ; js  = bxH%idxS(2) ; je  = bxH%idxE(2)
Isq = bxQ%idxS(1) ; Ieq = bxQ%idxE(1) ; Jsq = bxQ%idxS(2) ; Jeq = bxQ%idxE(2)
ksc = bxH%idxS(3) ; kec = bxH%idxE(3)      ! not ks/ke: see §9 #18
```

**Memory mode.** `G%IscB` is `isc` in non-symmetric and `isc-1` in symmetric memory, so growing an
h-box by one matches a B-index range in only one mode. Build staggered boxes from `G`'s B-indices
(`bxQ%set(idxS=[G%IscB,G%JscB,ks], idxE=[G%IecB,G%JecB,ke])`), never by growing.

**Tiles.** A root that loops over tiles allocates its boxes once and re-`set`s them per tile:

```fortran
call bxH%safe_alloc(ndims=3)
do k_start = 1, nz, nkblock
  k_end = min(k_start+nkblock-1, nz)
  call bxH%set(idxS=[G%isc,G%jsc,k_start], idxE=[G%iec,G%jec,k_end])
  call tile_body(..., bxH, ...)
enddo
call bxH%free()
```

In the tile body, scratch may take its k-bounds from the box (legal: a specification expression
may use a non-optional, non-`intent(out)` dummy, declared first) and be indexed by the real `k`:
`real, dimension(SZIB_(G),SZJB_(G),bxH%idxS(3):bxH%idxE(3)) :: q`. Full-k outputs written per
tile must be `intent(inout)`, or other tiles' layers are lost.

**Persistent grid containers `grid_core` / `grid_OBC`** (`src/core/MOM_grid_containers.F90`):
- Copies of `ocean_grid_type` arrays, built once in `initialize_MOM` and freed in `MOM_end`.
  `grid_core` holds every non-OBC array a planned tree reads; `grid_OBC` holds `OBCmaskCu
  IdxCu_OBCmask OBCmaskCv IdyCv_OBCmask`. Members are `RealArray_t` named as in `G`; each
  structure carries `C`, its `bind(C)` mirror (`grid_core_C`/`grid_OBC_C`), built under
  `#ifdef _TIM`.
- `Gcore => grid_core(G)` / `Gobc => grid_OBC(G)` return pointers; FATAL if not built or bounds
  don't match `G`.
- In a kernel: `type(grid_core_type), intent(in) :: Gcore`, a pointer per field, then
  `call Gcore%IdxCu%view(IdxCu)`; the body uses `IdxCu` for `G%IdxCu`.
- **Never** `%free()` a member, never `alloc(source=G%<f>)` per call.
- New members are appended in the same order in all five places (both types, init, `to_c`, end),
  comments copied word for word from `MOM_grid.F90`. Adding a field or group is the user's call.
- A tree root gets whole-column boxes (`bxH0`, `bxQ0`) from its cap and builds U/V tile boxes by
  combining their i/j ranges.

On the C++ side each array's own `lb` positions its `Array4` (protoMOMxx `generate_amrex_code`
lessons.md §7 #12); mixed staggerings break silently otherwise.

---

## 6. Intent-driven copy rules

| Original intent | Allocate | Copy in | Copy back (`copy2F`) | Container dummy intent |
|---|---|---|---|---|
| `intent(in)` | yes | yes — `source=x` | no | `intent(in)` |
| `intent(out)` | yes | **yes by default** — `source=x` | yes | `intent(inout)` |
| `intent(inout)` | yes | yes — `source=x` | yes | `intent(inout)` |

`intent(out)` means "undefined on entry", not "written in full". Grid-shaped arrays include the
halo; a `Box_t`-scoped write doesn't. Without `source=` the unwritten part is garbage that
`copy2F` copies over the caller's real data (§9 #2). Allocate-only (`alloc` without `source`, or
`allocView`) only when the callee provably writes the full `LBOUND:UBOUND`. In doubt, copy in.

### Optional dummies

A converted `optional` array dummy keeps `optional`. In the callee:

```fortran
  type(RealArray_t), optional, intent(in) :: hin_a  !< Initial layer thickness [H ~> m or kg m-2].
                                                    !! If hin is absent, h is also the initial thickness.
  real, dimension(:,:,:), contiguous, pointer :: hin

  nullify(hin)                                  ! §9 #16
  if (present(hin_a)) call hin_a%view(hin)      ! guard the view
  ...
  if (present(hin_a)) then                      ! every present(hin) renamed (§9 #8a)
```

Uses of the pointer stay inside the same guard. At the call site, if the caller's own source is
optional, allocate conditionally and branch the call — a container local is always present
(§9 #8b):

```fortran
  if (present(hin)) then
    call hin_a%alloc(lb=LBOUND(hin), ub=UBOUND(hin), source=hin)
    call continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a, hin_a=hin_a)
  else
    call continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a)
  endif
  call hin_a%free()      ! safe even if never allocated (§4.5)
```

Past two such optionals the branches explode: stop and ask (or see
`convert_optional_args_to_containers`). Optional scalars are never containerized.

---

## 7. Naming and comments

| Thing | Convention | Example |
|---|---|---|
| Subroutine | unchanged | `zonal_mass_flux` |
| Container dummy or call-site local | name + `_a` | `h_in` → `h_in_a` |
| `%view` pointer in the body | the original name | `h_in` |
| Grid field (persistent copy) | `%view` of `Gcore%<f>`, pointer named `<f>` | `G%mask2dT` → `mask2dT` |

The pointer keeps the original name so the math needs no edit.

**Doc comments.** Every existing `!<` comment survives verbatim (wording, units, punctuation),
including when a two-line declaration collapses to one — record it before rewriting the
declaration (see the STANDING RULE). Every new dummy gets a new `!<` comment describing its
physical meaning (a missing one breaks the Doxygen check). Comments not attached to a dummy stay
as they are.

**Line limit.** Keep every added or changed line within the file's established width and never
over 132 characters for code (MOM6 CI, `.testing/trailer.py`). Wrap with `&` or `!!`; never
shorten a comment, or drop or rename an argument, to fit. The same CI check rejects trailing
whitespace on any line; removing a statement from a `;`-joined line is the usual way to leave
some. Every skill's format check covers both (`grep -nE '[[:space:]]+$'` must find nothing).

**No commentary about the conversion.** Never add comments such as `! Containers for X` or
`! Must be nullified`; the diff is mechanism, not narration.

---

## 8. Derived-type fields and shrinking argument lists

Classify every `<dummy>%<field>` reference in the body (`convert_array_containers` Step 2 has the
full table):
- **`G` array** → a `%view` of the persistent `grid_core`/`grid_OBC` copy (§5).
- **Other array with a persistent container copy** (a bundle inside a control structure) → a
  `%view` of it. An array with no persistent copy is a question for the user; never a per-call
  `alloc(source=…)` dummy for data that has one.
- **Index range/extent** → bound of a `Box_t` dummy.
- **Scalar** (`GV%Angstrom_H`, `CS%upwind_1st`) → plain scalar dummy, derived at the call site.
- **Needed whole by a still-raw callee** → keep the derived type; don't force it.

Drop a derived-type dummy only when nothing in the body references or forwards it (§1).

---

## 9. Recurring pitfalls

1. **Don't rename the subroutine.** No `_fortran`, no wrapper (§10).
2. **Missing `source=` on a "pure output" corrupts its halo** (§6). Real CI failure:
   `zonal_mass_flux`'s `uh`.
3. **`grow`/`shrink`/`growLo`/`growHi`/`growBy` return a new box you must `free`** — else a leak
   per call, inside timestep loops.
4. **`set` before `safe_alloc` is FATAL.**
5. **`copy2Array` does not allocate** (§4.4).
6. **Written container dummies are `intent(inout)`, never `intent(out)`** — `intent(out)` resets
   the pointer components and discards the caller's allocation.
7. **Container extents must match the array's.** Copies loop over `this%shape`; a mismatch reads
   or writes out of bounds instead of failing.
8. **`view` rank mismatch FATALs at run time**, not compile time.
8a. **Rename every `present(<orig>)` to `present(<orig>_a)`, and guard the `%view`.** A missed
    rename fails to compile; a missed guard fails at run time.
8b. **Never pass an unallocated container to an `optional` dummy** — it reads as present. Branch
    the call when the caller's source is optional (§6).
9. **Update every call site in the same change** — `grep -rn "call[[:space:]]\+<name>(" src/
   config_src/`, plus continuation lines; check each position's container/raw kind, not just
   count and order.
10. **`elemental` callees keep their scalar signatures** (`flux_elem`, `flux_elem_OBC`); pass them
    `%view` pointer elements.
11. **A derived type with raw array members is not a container** (`BT_cont_type`). Converting a
    routine that passes it doesn't convert its members; leave it unless asked
    (`create_shadow_container_type`).
12. **Array `%to_c` is TIM-only** (FMS2 has it commented out): every array `%to_c()` sits inside
    `#ifdef _TIM`; `Box_t%to_c` exists in both. A conversion emits no `to_c`.
13. **`IntArray_C` is not public** in either layer; a caller can't declare `to_c_Int`'s result.
14. **Default `real` is assumed `real64`** (`copy2*` takes default `real`, `alloc*` takes
    `real(kind=real64)`). Never mix kinds across this.
15. **Free every container on every path**, early returns included.
16. **`nullify` every conditionally-`%view`'d pointer local before the first guard.** Its
    association status on entry is undefined, not disassociated; skipping this caused silent
    corruption (`present_uhbt_or_set_BT_cont`). Never `=> null()` in the declaration (implies
    `SAVE`). A nullified guarded pointer may be forwarded unconditionally to a still-raw callee's
    plain (non-pointer, non-allocatable) `optional` array dummy, where it reads as absent —
    `convert_optional_args_to_containers`.
17. **Never reference a container (`%associated()` or any component) inside an
    `!$omp target`/`teams`/`parallel`/`loop` construct.** Compute a plain scalar before it
    (`use_uhbt = uhbt_a%associated()`) and use only that. The implicit capture caused an ASan
    `bad-free` in `zonal_flux_adjust`.
18. **Names are case-insensitive.** A new `ke` clashes with an array `KE`; gfortran reports
    `Symbol 'ke' ... already has basic type of REAL` plus cascade errors elsewhere. `grep -niw` each
    new name before adding it, check per subroutine for case-insensitive duplicates, and read the
    first error of a failed build. Use `ksc`/`kec`.

---

## 10. Verification and the other skills

Verification (procedure: `convert_array_containers` Step 7): the diff is declarations, `%view`
calls and call-site blocks, with nothing inside the math; every call site matches; the code uses
only the shared API (§4), so it must build under both FMS2 and TIM; results are bit-identical.
Building is the user's job unless they allow it in this conversation (`~/.claude/CLAUDE.md`);
never call a conversion verified on inspection alone.

`generate_cpp_bridge` layers on afterwards: `*_fortran` rename, dispatcher, `getenv_mode`,
capture mode, `#ifdef _TIM` arm, `bind(C)` interface; `generate_amrex_code` writes the C++.

A raw-array subroutine with a `<name>_fortran` sibling is a bridge shim: convert it like any
other, but never edit between `select case (mode)` and `end select` (except renaming a reference
to a dropped derived type) or any `bind(C) interface … end interface`. Marshalling around them is
ordinary call-site work (§3). Needing more means the task has drifted into bridge territory.
