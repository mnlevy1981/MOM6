---
name: create_shadow_container_type
version: "0.3.2"
description: Build a container-based shadow of a shared MOM6 derived type (one whose fields are plain allocatable real arrays, but which is also used directly by code outside the call tree being bridged) so that one call tree can become container/bridge-ready without converting the shared type itself. Worked example in this codebase: BT_cont_type is used both by the continuity() call tree (MOM_continuity_PPM.F90) and directly by MOM_barotropic.F90; BT_cont_container_type shadows it inside continuity() only, leaving MOM_barotropic.F90 and MOM_variables.F90's BT_cont_type completely untouched. Use this when convert_array_containers doesn't apply because the type in question is shared outside the subroutine/call-tree being converted, so converting it wholesale would ripple into unrelated code you don't want to touch.
user-invocable: true
argument-hint: <work-directory> <shared-type-name> <call-tree-entry-point> [--enable_git_commit] [--disable_git_commit]
---

# Build a container-based shadow of a shared derived type

## Why this exists

`convert_array_containers` converts a subroutine's own dummy arguments.
That works cleanly when the subroutine "owns" the data it's converting.
It breaks down for a **shared derived type**: a type whose fields are
plain `real, allocatable` arrays (not yet containers), passed as a
single dummy of that type, but which is also read or written directly
by code entirely outside the call tree you're trying to make
container/bridge-ready.

**Worked example in this codebase.** `BT_cont_type`
(`src/core/MOM_variables.F90`) has 14 real allocatable fields plus 2
`group_pass_type` halo-pass fields. It is a dummy of `continuity()`
(`MOM_continuity_PPM.F90`'s call tree — the one being bridged) *and* is
accessed directly at 53 call sites in `MOM_barotropic.F90`, and passed
whole into `btcalc`/`btstep` from `MOM_dynamics_split_RK2[.b].F90`.
Converting `BT_cont_type` itself to containers would force every one of
those 53 accesses, plus the RK2/RK2b call sites, to change too — none of
which has anything to do with the bridging work at hand.

**The fix:** build a separate, container-based *shadow* type, local to
the call tree that needs it. Construct it from the real struct at the
call tree's own entry point, forward the shadow through every
descendant in place of the real type, and copy the results back into
the real struct before the entry point returns. The real shared type,
every *other* subsystem that touches it, and every subroutine outside
the call tree are never touched.

**`ocean_grid_type` is already handled differently:** it is static after init, so its arrays
have *persistent* copies (`grid_core`/`grid_OBC` in `src/core/MOM_grid_containers.F90`, see
`array_container_lessons` §5), built once rather than shadowed and copied back per call. Don't
shadow `G`.

This is a **pilot pattern**, not a one-off — the same technique
generalizes to any shared type too big or too widely used to convert
wholesale (e.g. `ocean_grid_type`/`verticalGrid_type`). Read the Steps
below as a checklist for any such type, not just `BT_cont_type`.

This skill assumes `array_container_lessons` and (if a bridge already
exists somewhere in the call tree) `cpp_bridge_lessons` have already
been invoked this session, for the container API (`%alloc`/`%view`/
`%copy2F`/`%free`/`%associated()` signatures) and the `%to_c()`/
`bind(C)` conventions. Neither is reproduced here.

**Not `create_config_bundle_type`.** If the type in question is
*private* to the call tree you're working on (no use outside it), you
don't need a shadow-and-copy-back dance at all — restructure the real
type's fields in place instead. Use that sibling skill there; use this
one only when the type is genuinely shared.

## Precondition

The shared type's fields that matter to this call tree must be plain
`real`/`integer`/`logical`, allocatable or automatic arrays — the same
element types `convert_array_containers` already knows how to
containerize (`RealArray_t`/`IntArray_t`/`LogicalArray_t`). If a field
is itself a derived type (e.g. `BT_cont_type`'s `group_pass_type` halo
fields), it is out of scope for the shadow — see Step 1's "drop
unused fields" rule below.

Before starting, read every field of the real type and grep the call
tree for `<dummy>%<field>` to find which fields the call tree actually
touches. A field the call tree never reads or writes does not belong in
the shadow at all, regardless of whether it's a plain array or a nested
derived type.

## Steps

### 1. Define the shadow type (self-contained, zero risk, independently committable)

In the module that hosts the call tree's leaf kernels (not the module
that defines the real shared type — the shadow is scoped to *this* call
tree, and lives with the code that uses it), add:

```fortran
!> A container-based shadow of <RealType>, used only within the <entry_point>() call tree so
!! that its <field-group> fields can be bridge-ready without converting the shared <RealType>
!! itself (which is also used directly by <other consumer(s)>).
type, public :: <RealType>_container_type ; private
  type(RealArray_t) :: field1_a, field2_a, ... ! one per field the call tree touches
contains
  procedure, public :: build_from  => <RealType>_container_build_from
  procedure, public :: copy_back   => <RealType>_container_copy_back
  procedure, public :: associated  => <RealType>_container_associated
end type <RealType>_container_type
```

- **`build_from(real_struct)`** — alloc + copy-in each field from the
  real struct, exactly like a `convert_array_containers` call-site
  marshalling block (lessons §6), but written once as a type-bound
  procedure instead of repeated at every call site.
  - **Respect the real type's own allocated-together grouping.** If the
    real type allocates several fields together in one init-time
    routine (check that routine, don't guess), guard that whole group
    with a single `allocated(real_struct%<one_field_from_the_group>)`
    check, mirroring the real type's own invariant instead of inventing
    a per-field one. `BT_cont_type`'s 12 `FA_*`/`*BT_*` fields are
    always allocated together (per `alloc_BT_cont_type`) and need no
    guard; `h_u`/`h_v` are allocated together only if the caller
    requested `alloc_faces` at init time, and are left unassociated in
    the shadow (and therefore null-safe downstream) when absent.
- **`copy_back(real_struct)`** — the reverse: `%copy2F` each field into
  the real struct, then `%free()` it. Guard exactly the same field
  groups the same way.
- **Fields allocated independently.** If the real type allocates each field on its own,
  so that there is no allocated-together group (CorAdCalc's `CorAdv_AD_type` shadow of
  `accel_diag_ptrs`: `gradKEu`/`gradKEv`/`rv_x_u`/`rv_x_v`, allocated by several modules):
  - guard each field separately in `build_from` (`associated(AD%x)`) and `copy_back`
    (`x_a%associated()`);
  - leaves test and `%view` each field under its own `%associated()`, because `%view` on an
    unallocated container is FATAL;
  - omit the whole-shadow `associated()` below; there is no valid sentinel, and none should be
    invented.
- **`associated()`** — a `pure function` returning whether this shadow
  was ever built. Pick **one field from the group that's always
  allocated together** as the presence sentinel (do not add a separate
  boolean flag) — `BT_cont_container_type` uses
  `this%FA_u_EE_a%associated()` since the 12 always-together fields
  make any one of them equally valid as the signal.

Drop any field the call tree never touches (per the precondition check)
— do not carry a derived-type field like `group_pass_type` into the
shadow just because the real type has it; if the call tree needs a halo
pass, that stays on the real struct, done by the code that already owns
it.

This stage touches no existing subroutine — it's pure dead-code
addition and can be reviewed/committed on its own before Step 2 begins.

### 2. Thread the shadow through the call tree (one atomic change)

The shadow's type has to match at every call boundary through the tree,
so this can't be split level-by-level — an intermediate state where one
level passes the shadow into a level still expecting the real type is a
straight compile error. All of the following land together, walking the
tree from the entry point down:

1. **The call tree's entry point** (the outermost wrapper — its own
   public signature does not change; external callers keep passing the
   real struct). In its body, alongside the wrapper's other optional
   -argument marshalling: build the shadow via `%build_from()`,
   guarded exactly like the real dummy's own optionality (`if
   (present(real_dummy)) then ; if (associated(real_dummy)) call
   shadow%build_from(real_dummy) ; endif`). Forward the shadow into the
   next level in place of the real dummy. After that call returns,
   mirror the wrapper's existing copy-back pattern: `if
   (shadow%associated()) call shadow%copy_back(real_dummy)`.

2. **Every intermediate level** — retype its own dummy from
   `type(RealType), optional, pointer` (or however the real type
   arrived there) to `type(RealType_container_type), intent(inout) ::
   shadow_a`, renamed with the container-name suffix already
   established in this codebase (`_a`). This is a pure rename/retype at
   levels that only forward the struct without touching a field — no
   logic change. **Every internal "argument is absent" sentinel** built
   at a level that calls a leaf without one (e.g. a mode that never
   uses this feature) switches from simply omitting the optional
   keyword to passing a never-`%build_from`'d shadow instance — its
   `%associated()` reads false, which is exactly the "absent" signal the
   leaf now expects.

3. **A level that already has a `bind(C)` bridge and only needs the
   struct as an opaque pointer** (never touches a field directly) —
   confirm this really is opaque (grep the kernel's own body for
   `<dummy>%<field>`; if it's genuinely never dereferenced there, its
   `bind(C)` interface already declares it `type(c_ptr), value`, built
   via `associated(dummy)` → `c_loc(dummy)`/`c_null_ptr`.
   **`c_loc`/`associated` don't care what Fortran type the pointer
   target is** — retyping the dummy from the real struct to the shadow
   at this level needs **zero changes to the `bind(C)` interface or the
   AMReX arm**, only the variable name in the `c_loc`/`associated` call.
   Confirm this by re-reading the AMReX arm after the rename: it should
   be textually identical apart from the name. If a `target` attribute
   is needed on the shadow dummy for `c_loc()` to be legal there (it
   wasn't on the real struct, since that arrived as a `pointer`, not a
   plain `intent(inout)` value), add it — this is the one attribute
   change this case does require.

4. **A leaf that writes real-struct fields directly** — replace every
   `real_struct%field(i,j) = ...` with a plain local name obtained via
   `call shadow%field_a%view(field)` once near the top of the body,
   then leave every element-wise expression byte-identical (rename the
   left-hand side only, exactly like an ordinary
   `convert_array_containers` dummy conversion). Guard the `%view`
   the same way the leaf's own "is this feature active" flag already
   guards the real-struct writes today — do not add a new guard that
   wasn't there before.

5. **A leaf that currently round-trips one field through a throwaway
   local container purely to hand it to an already-converted
   descendant** — this collapses for free once the field lives in the
   shadow: **Fortran matches actual arguments by type/rank/shape, not
   variable name**, so `shadow%field_a` can be passed directly into
   the descendant's dummy slot. The old `%alloc(source=real_struct
   %field)` / call / `%copy2F(real_struct%field)` / `%free()`
   round-trip disappears entirely — this is the same elimination
   `hoist_container_marshalling` §3 describes for two separately-built
   containers holding the same value, except here it falls out
   automatically because the shadow's own field *is* that value, with
   no separate step needed.

Verify after this step, before moving on: grep every leaf and
intermediate subroutine touched for a remaining `<realtype>%<field>` or
bare real-typed dummy reference — zero should remain outside the entry
point's own `build_from`/`copy_back` calls.

### 3. Add a `bind(C)` mirror, only if and when a leaf needs structured field access at the C boundary

Skip this step entirely if every leaf that crosses into C treats the
struct as opaque (Step 2 item 3) — do not add a mirror type nobody
calls.

When a leaf's bridge genuinely needs the real fields on the other side
(e.g. `set_zonal_BT_cont`'s C++ counterpart needs to write into
specific face-area fields, not just receive an opaque handle):

- **Define the mirror next to the real shared type, in the module that
  defines it — never in the shadow's own module or a generic
  container-primitives module.** Placing it in the shadow's module
  produces a real circular-`use` compile error (full mechanism in
  Pitfalls, below) — and it's the module a future consumer of the real
  type would look for it in anyway.
  ```fortran
  !> bind(C) mirror of <RealType>'s N real fields, one RealArray_C per field (same order as
  !! <RealType> itself). Lives here, next to <RealType>, rather than in <shadow-module> or a
  !! container-primitives module -- it mirrors this real, allocatable-field type at the bridge
  !! boundary, not a container type.
  type, bind(C) :: <RealType>_C
    type(RealArray_C) :: field1 !< See <RealType>%field1.
    ...
  end type <RealType>_C
  ```
  Export it `public` from that module.
- Every `bind(C)` interface elsewhere that needs the mirror type
  imports it with `use <RealType's module>, only : <RealType>_C` inside
  its own interface body — never `use <shadow's own module>`, even if
  that's where the interface block physically sits.
- Add a **standalone function**, not a type-bound procedure, that
  builds the mirror from the shadow by calling each field's own
  `%to_c()`:
  ```fortran
  #ifdef _TIM
  function <RealType>_container_to_c(shadow_a) result(cdesc)
    type(<RealType>_container_type), intent(in) :: shadow_a
    type(<RealType>_C) :: cdesc
    cdesc%field1 = shadow_a%field1_a%to_c()
    ...
  end function <RealType>_container_to_c
  #endif
  ```
  It must be a free function, outside the shadow type's own
  unconditionally-compiled `contains` block, because `%to_c()` itself
  only exists under the TIM infra layer (FMS2 has it fully commented
  out) — every other `%to_c()` call in the file already lives inside
  an `#ifdef _TIM` arm for the same reason (`cpp_bridge_lessons` §3.1,
  §10). **Put the doc comment for this function *before* the `#ifdef
  _TIM` line, not between it and the `function` line** — see Pitfalls
  below for why the more natural-looking placement produces a doxygen
  false positive that is easy to "fix" incorrectly once.

### 4. Add capture-mode support to the shadow type

If any leaf's shim records inputs/outputs via `io_recorder` for offline
validation, the shadow type needs to be recordable too, or a
capture/replay run silently drops it.

- Add `%write_binary()`/`%read_binary()` **directly to the shadow
  type**, delegating to each field's own `%write_binary()`/
  `%read_binary()` in one fixed order (the same order used everywhere
  else in the type). These are **unconditionally compiled — no
  `#ifdef _TIM` needed** — they never touch `%to_c()` or the AMReX
  shadow buffer, only `%data` via each field's own binary I/O, which is
  available identically under both infra layers.
- **`io_recorder` (the framework's capture module) cannot gain a
  generic `add`/`get` overload for the shadow type directly** — that
  would require the framework module to `use` the module the shadow
  type lives in, which already `use`s the framework module: a circular
  dependency. Instead, add two **free subroutines in the shadow type's
  own module**, `record_<name>`/`restore_<name>`, that compose the
  framework module's own already-public primitives directly (its
  `add_entry`/`find_entry`/`entries`/`unit_bin` — check they're public
  by confirming the framework module has no blanket `private`
  statement covering them; mirror exactly what its own
  `add_realarray`/`get_realarray` do internally):
  ```fortran
  subroutine record_<name>(rec, name, shadow_a)
    type(io_recorder),               intent(inout) :: rec
    character(*),                    intent(in)    :: name
    type(<RealType>_container_type), intent(in)    :: shadow_a
    integer(kind=int64) :: pos
    inquire(unit=rec%unit_bin, pos=pos)
    call rec%add_entry(name, '<RealType>_container_type', pos)
    call shadow_a%write_binary(rec%unit_bin)
  end subroutine record_<name>

  subroutine restore_<name>(rec, name, shadow_a)
    type(io_recorder),               intent(inout) :: rec
    character(*),                    intent(in)    :: name
    type(<RealType>_container_type), intent(inout) :: shadow_a
    integer :: idx
    integer(kind=int64) :: pos
    idx = rec%find_entry(name)
    if (idx < 0) call MOM_error(FATAL, "restore_<name>: variable "//trim(name)//" not found")
    if (trim(rec%entries(idx)%type_name) /= '<RealType>_container_type') &
      call MOM_error(FATAL, "restore_<name>: variable "//trim(name)//" type mismatch")
    pos = rec%entries(idx)%offset
    read(rec%unit_bin, pos=pos)
    call shadow_a%read_binary(rec%unit_bin)
  end subroutine restore_<name>
  ```
  This needs a `use, intrinsic :: iso_fortran_env, only : int64` at the
  host module's top level if not already present.
- Call `record_<name>`/`restore_<name>` from every leaf shim's capture
  arm that touches the shadow, using the same `_before`/`_after` naming
  convention as every other capture entry (`cpp_bridge_lessons` §13).
- If an external test harness needs to call these directly (rather than
  only from inside a shim's capture arm), add them to the host module's
  `public ::` list — this is a deliberate, separate decision from
  making the shadow type itself public (it already is, from Step 1),
  since a free subroutine needs its own explicit export.

### 5. Verify

- **Everything outside the call tree is untouched.** Diff the module
  that defines the real shared type, and every other subsystem that
  uses it directly (in the worked example: `MOM_variables.F90`,
  `MOM_barotropic.F90`, the RK2/RK2b dynamics files) against HEAD —
  must be empty except for the Step 3 mirror addition (which lands in
  the real type's module, additively, and changes nothing else there).
- **The entry point's own public signature is unchanged**, and every
  one of its external callers' call sites is byte-identical — the
  shadow is built and torn down entirely inside the entry point's own
  body.
- **Re-read every already-bridged leaf's `bind(C)` interface and AMReX
  arm after the rename** (Step 2 item 3) to confirm it is textually
  identical apart from the dummy's name — no interface signature or
  bridge-logic change should have been needed there.
- Grep every touched leaf/intermediate subroutine for a lingering
  `<realtype>%<field>` or a raw dummy of the real type — none should
  remain outside the entry point.
- Confirm the shadow's `build_from`/`copy_back` preserve every
  allocated-together-or-not-at-all grouping from the real type's own
  init routine (Step 1).
- Doc-comment count (`!<`) and line length (≤100) checks on every added
  or renamed line, same as `convert_array_containers` Step 9.
- If a `bind(C)` mirror was added (Step 3): confirm its `#ifdef _TIM`/
  `#endif` pair has no doc comment sandwiched between the `#ifdef` and
  the function/subroutine it guards (Pitfalls below); run a depth-
  tracking check across the whole file (`#ifdef`/`#ifndef`/`#if` = +1,
  `#endif` = -1, must never go negative and must end at 0) rather than
  trusting a raw count of the two tokens.
- No Fortran compiler locally — report the build as unrun, same as
  every other skill in this family.

## Pitfalls (each one a real bug hit while building this pattern)

- **A module cannot `use` itself, even from a nested `interface` body.**
  Placing the `bind(C)` mirror type in the shadow's own module and then
  having a `bind(C)` interface in that *same* module `use` it by that
  module's name fails to compile ("Cannot open module file") — the
  `.mod` file doesn't exist yet mid-compilation. Fix: the mirror lives
  in the real shared type's module (Step 3), never the shadow's.
- **A type used only inside a `bind(C)` interface body isn't visible to
  the module's own contained shim procedures.** An `interface` block is
  a separate scoping unit — it does not get host association from the
  module's top-level `use` statements, and neither do the module's own
  `contains`-ed procedures get whatever the interface body imported.
  If a shim declares a local of the mirror type (e.g. to build it
  before a bridge call), that type must be imported at the **module's
  own top-level** `use` line too, not only inside the interface body —
  a real compile error ("derived type ... being used before it is
  defined") results otherwise.
- **A doc comment sandwiched between `#ifdef _TIM` and the code it
  guards produces a "More #endif's than #if's found" doxygen warning,
  even when the `#ifdef`/`#endif` count is textually balanced.** This
  is not caught by counting the two tokens; it took a custom
  depth-tracking script to isolate, because every *other* `#ifdef
  _TIM` pair in the same file had the guarded code immediately after
  the directive with no comment in between, and this was the only one
  that didn't. Fix: put the doc comment *before* `#ifdef _TIM`, matching
  every other guarded pair's shape exactly, rather than between the
  directive and the function signature it appears to be documenting.
- Rephrasing a doc comment to fix the sandwiched-comment warning above
  can accidentally fix a *different*, textually-similar-looking warning
  too ("explicit link request to `ifdef` could not be resolved") if
  the comment's prose literally contains the substring `#ifdef` —
  doxygen's scanner treats that text as an autolink attempt. Don't
  assume the two warnings share one root cause just because they
  disappear at the same time; verify each independently (a link-target
  warning is a wording problem, the endif-count warning is a structural
  placement problem, and this pattern can trigger either, both, or
  neither depending on exactly how the comment is phrased).

## Versioning marker

Every Fortran file this skill creates or modifies gets a `!!SKILLS: 0.3.2`
marker line — the shared version number for this whole skill family,
not just this one skill (bump every skill file's `version:` field and
this marker in lockstep when any of them changes in a way that affects
generated code, never per-skill). If the file doesn't already have one,
add it as its own line immediately after the file's existing license/
header comment block, before the `module` statement. If it already has
one, update it to the current version rather than adding a second line.
Deliberately grep-able (`grep -rn "!!SKILLS:"`) and meant to be stripped
later, once these markers are no longer useful.

## Hard rules

- Never skip the `!!SKILLS: 0.3.2` marker on a file this skill touches,
  and never add a second marker line if one already exists — update it
  in place instead.
- Never convert the real shared type itself, and never touch any
  subroutine or module outside the call tree that owns the shadow —
  that is the entire reason this pattern exists instead of
  `convert_array_containers`.
- Never add a field to the shadow that the call tree doesn't actually
  read or write, even if the real type has it. Never carry a
  derived-type field (e.g. a halo-pass structure) into the shadow at
  all.
- Never invent a presence flag for the shadow type — always derive
  `%associated()` from an existing field that's part of the
  always-allocated-together group (Step 1).
- Never place a `bind(C)` mirror type in the shadow's own module, or in
  a generic container-primitives module — it belongs next to the real
  shared type it mirrors (Step 3), and this is not merely a style
  preference: the wrong location produces a real circular-`use` compile
  error under the conditions described in Pitfalls.
- Never skip Step 2's per-level check of whether a leaf's `bind(C)`
  interface treats the struct as opaque (item 3) versus needs
  structured field access (item 4/Step 3) — assuming the wrong one
  either produces unnecessary bridge-signature churn or misses a case
  where fields genuinely need to cross into C.
- Do not add `%to_c()` (or the standalone `_to_c()` function from Step
  3) unguarded by `#ifdef _TIM` — it depends on a type-bound procedure
  that does not exist under the FMS2 infra layer at all.
- Do not claim a build passed that was never run, and do not attempt to
  install a Fortran compiler.
- No commentary in the source about the migration itself (why a field
  became a container, why a guard is shaped the way it is) beyond what
  a genuinely new or changed doc comment needs to describe the
  physical quantity — this codebase documents ocean physics, not its
  own refactoring history, same standing rule as every other skill in
  this family.

## Commit gating

Same as `convert_array_containers`: `--enable_git_commit` /
`--disable_git_commit` override; otherwise follow
`~/.claude/preferences.json`'s `git_commit_and_push` key (`"auto"` runs
the commit step; `"manual"`, a missing file, a missing key, or
unparseable JSON all mean skip it and report modified files for manual
commit). Branch name:
`claude_<lowercased_realtype-name>_shadow_container`.

## Output to the user on success

1. The shared type shadowed, the call tree's entry point, and which of
   its fields made it into the shadow (and which were dropped as
   unused, and why).
2. Every level threaded through Step 2, tagged by which of the three
   cases applied (opaque-pointer rename, direct-field `%view`, or
   round-trip elimination).
3. Whether a `bind(C)` mirror and capture-mode support were added
   (Steps 3–4), and, if so, where the mirror type lives and which
   `record_*`/`restore_*` helpers were exported.
4. Confirmation (via diff) that the real shared type's own module and
   every other subsystem using it directly are untouched.
5. Build status — explicitly unrun if no compiler is available, per
   standing project policy.
6. Whether committed, or the list of modified files for manual commit.
