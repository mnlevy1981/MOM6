---
name: create_shadow_container_type
version: "0.4.0"
description: Build a container-based shadow of a shared MOM6 derived type -- one whose plain allocatable array fields are also used directly by code outside the call tree being converted -- so that call tree becomes container/bridge-ready without touching the shared type. The shadow is built from the real struct at the tree's entry point, threaded down in its place, and copied back before return (e.g. BT_cont_container_type shadows BT_cont_type inside continuity() only). Use when the type is shared outside the tree; for a private, module-local type use create_config_bundle_type instead.
user-invocable: true
argument-hint: <work-directory> <shared-type-name> <call-tree-entry-point> [--enable_git_commit] [--disable_git_commit]
---

# Build a container-based shadow of a shared derived type

`convert_array_containers` converts a subroutine's own dummies. A **shared** type (plain
`real, allocatable` fields, also read or written directly by code outside the tree) can't be
converted wholesale without rippling into unrelated code. Instead, build a container shadow
local to the tree: construct it from the real struct at the entry point, forward it in place
of the real type, copy back before return. The real type and everything outside the tree stay
untouched. Example: `BT_cont_type` (`MOM_variables.F90`) is a dummy of `continuity()` and is
also used at 53 sites in `MOM_barotropic.F90`; `BT_cont_container_type` shadows it inside
`continuity()` only.

- **Don't shadow `G`.** `ocean_grid_type` is static after init and has persistent copies
  (`grid_core`/`grid_OBC`, `array_container_lessons` §5).
- **Private type?** Use `create_config_bundle_type` and restructure it in place; no shadow.
- Container API and `%to_c()`/`bind(C)` conventions come from `array_container_lessons` and
  `cpp_bridge_lessons` (invoke once per session); they are not repeated here.

**Precondition.** Read every field of the real type and grep the tree for `<dummy>%<field>`.
Only plain `real`/`integer`/`logical` array fields the tree actually touches go into the shadow.
An untouched field, or a derived-type field (e.g. a `group_pass_type` halo pass), stays on the
real struct with the code that owns it.

## Steps

### 1. Define the shadow type (pure addition, committable on its own)

In the module hosting the tree's leaf kernels (not the real type's module):

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

- **`build_from(real_struct)`**: `%alloc(source=...)` each field, as in a call-site marshalling
  block (lessons §6). **`copy_back(real_struct)`**: `%copy2F` each field back, then `%free()`.
- **Guards mirror the real type's own allocation invariant**, read from its init routine, never
  invented. Fields allocated together share one `allocated(real_struct%<one field>)` guard (or
  none, if always allocated: `BT_cont_type`'s 12 `FA_*`/`*BT_*`; `h_u`/`h_v` only with
  `alloc_faces`, left unassociated otherwise). `copy_back` guards the same groups.
- **`associated()`**: a `pure function` returning whether the shadow was built, read from one
  field of an always-together group (`this%FA_u_EE_a%associated()`). Never add a boolean flag.
- **Independently allocated fields** (e.g. `CorAdv_AD_type` shadowing `accel_diag_ptrs`):
  guard each field on its own (`associated(AD%x)` / `x_a%associated()`), and have leaves `%view`
  each under its own `%associated()` (`%view` on an unallocated container is FATAL). Omit the
  whole-shadow `associated()`; there is no valid sentinel.

### 2. Thread the shadow through the call tree (one atomic change)

Every call boundary must agree, so all levels change together:

1. **Entry point**: public signature unchanged. Build the shadow under the real dummy's own
   optionality (`if (present(x)) then ; if (associated(x)) call x_a%build_from(x) ; endif`),
   pass it down, and afterwards `if (x_a%associated()) call x_a%copy_back(x)`.
2. **Intermediate levels**: retype the dummy to `type(<RealType>_container_type),
   intent(inout) :: <name>_a`; pure rename where it is only forwarded. Where a level used to
   omit the optional argument to signal "absent", pass a never-built shadow instead (its
   `%associated()` is false).
3. **A bridged level that treats the struct as opaque** (grep confirms no field access): its
   `bind(C)` interface already takes `type(c_ptr), value` via `associated()`/`c_loc()`, which
   don't care about the target type. Only the variable name changes, plus a `target` attribute
   on the dummy if `c_loc` now needs it; the interface and AMReX arm must stay otherwise
   textually identical.
4. **A leaf that writes fields directly**: `call <name>_a%field_a%view(field)` once near the top,
   under the leaf's existing feature guard (no new guard); element-wise expressions stay
   byte-identical.
5. **A leaf that round-trips a field through a throwaway container** for a converted callee:
   pass `<name>_a%field_a` directly (actuals match by type/rank, not name); the
   alloc/call/copy2F/free disappears (cf. `hoist_container_marshalling` §3).

### 3. Add a `bind(C)` mirror, only if a leaf needs field access at the C boundary

Skip it if every bridged leaf treats the struct as opaque. Otherwise:

- Define the mirror **in the real type's module**, next to the real type, and export it
  `public`. Never in the shadow's module (circular `use`, see Pitfalls) or a generic
  container-primitives module.
  ```fortran
  !> bind(C) mirror of <RealType>'s N real fields, one RealArray_C per field (same order as
  !! <RealType> itself).
  type, bind(C) :: <RealType>_C
    type(RealArray_C) :: field1 !< See <RealType>%field1.
    ...
  end type <RealType>_C
  ```
- `bind(C)` interfaces import it with `use <real type's module>, only : <RealType>_C`.
- The converter is a **free function**, not type-bound, guarded by `#ifdef _TIM` because
  `%to_c()` exists only under TIM (`cpp_bridge_lessons` §3.1). Its doc comment goes **before**
  the `#ifdef`:
  ```fortran
  !> Builds the bind(C) mirror of a <RealType>_container_type.
  #ifdef _TIM
  function <RealType>_container_to_c(shadow_a) result(cdesc)
    type(<RealType>_container_type), intent(in) :: shadow_a
    type(<RealType>_C) :: cdesc
    cdesc%field1 = shadow_a%field1_a%to_c()
    ...
  end function <RealType>_container_to_c
  #endif
  ```

### 4. Add capture-mode support to the shadow type

Needed if any leaf shim records via `io_recorder`; otherwise capture silently drops the shadow.

- Type-bound `%write_binary(unit)`/`%read_binary(unit)` delegating to each field's own, in one
  fixed order. Unconditionally compiled (no `%to_c()`, no `#ifdef`).
- `io_recorder` can't get a generic overload for the shadow (circular `use`). Add free
  subroutines in the shadow's module that compose its public primitives (`add_entry`,
  `find_entry`, `entries`, `unit_bin`; mirror `add_realarray`/`get_realarray`). Needs
  `use, intrinsic :: iso_fortran_env, only : int64`.
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
- Call them from every capture arm that touches the shadow, with the `_before`/`_after` naming
  (`cpp_bridge_lessons` §13). Add them to `public ::` only if an external harness needs them.

### 5. Verify, by script

- Diff of the real type's module and every outside consumer is empty, except the additive Step 3
  mirror. The entry point's signature and its external call sites are unchanged.
- No `<real dummy>%<field>` or real-typed dummy remains in the tree outside the entry point's
  `build_from`/`copy_back`.
- `build_from`/`copy_back` guards match the real type's init-time allocation groups.
- Every opaque bridged leaf's `bind(C)` interface and AMReX arm differ only by the dummy's name.
- Comments: every existing `!<` survives verbatim; line limit (lessons §7); `#if`/`#endif` depth-tracked across the file (never negative, ends at 0),
  and no doc comment between an `#ifdef _TIM` and the code it guards.
- **Build:** don't compile, build or run tests unless the user explicitly allows it in this
  conversation (`~/.claude/CLAUDE.md`). Report it as unrun, for both infra layers.

## Pitfalls (real compile/doxygen errors)

- **Circular `use`:** a mirror placed in the shadow's module and `use`d by name from a `bind(C)`
  interface in that same module fails ("Cannot open module file"). Put it in the real type's
  module.
- **Interface bodies are separate scopes:** a type imported only inside an `interface` body is
  invisible to the module's procedures. A shim with a local of the mirror type needs it on the
  module's top-level `use` too ("derived type ... being used before it is defined").
- **Doc comment between `#ifdef _TIM` and the function** gives doxygen "More #endif's than #if's
  found" even when balanced. Put it before the `#ifdef`.
- Comment prose containing the literal `#ifdef` gives a separate doxygen warning ("explicit link
  request to `ifdef` could not be resolved"). It is a wording problem, distinct from the
  placement one; verify each independently.

## Invariants (hard rules)

- Never convert the real shared type or touch anything outside the tree (Step 3's mirror is the
  one additive exception).
- Shadow holds only plain array fields the tree touches; never a derived-type field.
- `%associated()` comes from an existing always-together field, never a new flag.
- Decide opaque vs field-access per bridged level (Step 2.3 vs Step 3) by grep, not assumption.
- `%to_c()` and the `_to_c` converter only under `#ifdef _TIM`.
- Never reword, shorten or "fix" an existing comment; no comments about the refactor itself.
- Never claim an unrun build passed; never install a compiler.
- Exactly one `!!SKILLS: 0.4.0` marker per touched Fortran file, on its own line after the
  license header, before `module`; update an existing one in place. The version is shared by
  the whole skill family and bumped in lockstep.

## Commit gating

As `convert_array_containers` Step 8: `--enable_git_commit`/`--disable_git_commit` override;
otherwise `~/.claude/preferences.json` `git_commit_and_push` (`"auto"` commits; anything else,
missing or unreadable, skips and lists the files). Branch:
`claude_<lowercased_realtype-name>_shadow_container`. Never force-push.

## Output to the user

1. Type shadowed, entry point, fields included and fields dropped (and why).
2. Each level threaded and which Step 2 case applied.
3. Whether a mirror and capture support were added, where, and which helpers are public.
4. Step 5 results (including the untouched-outside diff) and build status; committed, or the
   files to commit.
