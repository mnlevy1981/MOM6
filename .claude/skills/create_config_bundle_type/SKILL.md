---
name: create_config_bundle_type
version: "0.4.0"
description: Cluster a private (module-local) MOM6 control-structure/config type's fields -- scalars and already-containerized arrays -- into small bind(C)-ready bundle types nested back into the original type, based on which fields actually travel together across the call tree (e.g. continuity_PPM_CS -> reconstruction_CS + transport_adjust_CS). Also gives init-time constant arrays persistent container copies (grid_core pattern) when the plan records it. Use instead of leaving individually flattened fields as the end state, when the type is private to the tree; for a type shared outside the tree use create_shadow_container_type instead.
user-invocable: true
argument-hint: <work-directory> <type-name> [<call-tree-entry-point>] [--enable_git_commit] [--disable_git_commit]
---

# Cluster a private control structure's fields into bind(C)-ready bundles

A control structure that can never be `bind(C)` (e.g. it holds a `type(diag_ctrl), pointer`)
gets flattened into individual dummies by `convert_array_containers` Step 2 wherever a bridged
leaf needs its fields. Left that way, the same few fields reappear in every signature and bridge
interface down the tree. Instead, trace which fields travel together and give each such cluster
its own small bundle type (plain fields, a `bind(C)` mirror, a converter), nested back into the
original type. The un-bindable fields stay on the original type, untouched.

- **Shared type?** If code outside the tree accesses its fields, stop and use
  `create_shadow_container_type`.
- Works on an already-flattened or an unflattened type; Step 6's per-leaf decision applies
  either way.
- Container API and `%to_c()`/`bind(C)` conventions come from `array_container_lessons` and
  `cpp_bridge_lessons` (invoke once per session).

## Field eligibility

- **Scalar** (any intrinsic type): eligible.
- **Array**: eligible only if already a container. Never convert a raw array here; run
  `convert_array_containers` first. Exception: init-time constant arrays, when the plan
  records it (Step 1b).
- **Pointer to a derived type, or anything never `bind(C)`-representable**: never bundled;
  stays on the original type.
- **Diagnostic handles** (`id_*`): normally excluded by the plan (they never cross to C++). If one
  also **gates computation** in a kernel (e.g. `if (CS%id_CAuS > 0 .or. CS%id_CAvS > 0)`), either
  report the leaf as needing the whole struct (so Phase 3 passes the gate as an argument), or add
  a logical member to the bundle (`do_CAS_diag`), set once at the end of init after diagnostics
  are registered, from the same test. The `id_*` handles stay on the struct for `post_data`.

A cluster may mix scalars and containers; the converter treats each member by its kind (Step 4).

## Precondition: the type is private to the tree

Grep the whole repo for `use <module>, only : <TypeName>` and `<module>%<TypeName>`. The type is
disqualified only if code outside the tree **accesses its fields**. Outside code that merely
declares and passes it is fine when the components are `private` (`type, public :: X ; private`)
— e.g. `CoriolisAdv_CS`, passed opaquely by the dynamics drivers. `continuity_PPM_CS` qualifies;
`BT_cont_type` (fields used by `MOM_barotropic.F90`) does not.

## Settle these decisions

Take the user's or the plan's values where given; otherwise:

1. **Bundle uniformly.** A subroutine needing only one or two fields of a cluster still takes the
   whole bundle. One exception: the hottest inner-loop leaf (`elemental`, per cell) may get a
   bare scalar unpacked by its caller at that one call site (`flux_elem`'s `vol_CFL`). Name it in
   the report; never apply it more widely.
2. **Clusters come from traced usage** (Step 2), never from field names or declaration order.
3. **No bundle for a field without a recurring companion.** A field used at one or two call
   sites, never with another field, stays a bare dummy (as `GV`/`US` fields in continuity did).
4. **Naming** (precedent: `MOM_continuity_PPM.F90`): type `<cluster>_CS`, mirror
   `<cluster>_CS_C`, converter `<cluster>_CS_to_c`; nested as a member named like its type
   (`type(reconstruction_CS) :: reconstruction_CS`); a leaf taking the bundle names the dummy
   `CS`, so its `CS%<field>` references need no edit. Exception: Step 1b's dummy.
5. **Stage multiple clusters**: each cluster is one atomic pass, since every call boundary for
   its fields must agree at once.

## Steps

### 0. Validate inputs

`<work-directory>` and `<type-name>` are required; `<call-tree-entry-point>` defaults to every
public subroutine in the type's module. Empty or `help`/`--help`/`-h`: print usage (from the
frontmatter) and stop. A missing work directory stops. Locate the declaration
(`grep -irn "^[[:space:]]*type.*::[[:space:]]*<type-name>\b"`): none is an error; more than one,
ask.

### 1. Classify every field

Read the full declaration and classify each field as above. A raw array you would cluster (and
that is not a Step 1b case) stops the run; report it for `convert_array_containers`.

### 1b. Persistent containers for init-time constant arrays (only when the plan records it)

For raw `allocatable` (or `ALLOCABLE_`) fields written only during init and only read by the
tree: metric products and coefficient fields (`hor_visc_CS`'s `dx2h`, `Kh_bg_xx`,
`Laplac2_const_xx`, …). Without this they would have to be copied into fresh containers at every
call, every timestep. User decision 2026-10-23 (`horizontal_viscosity` D4).

- **Pattern `grid_core`** (`src/core/MOM_grid_containers.F90`): keep every raw field as it is
  (init, `*_end` and Fortran-only code still use it) and add a container-only bundle nested in
  the type (e.g. `hv_metrics_CS`) holding one persistent copy per field. Members keep the raw
  field's name and its doc comment verbatim.
- **Build once at the very end of init**, after every helper that writes the arrays:
  `call CS%<bundle>%<f>%alloc(lb=LBOUND(CS%<f>), ub=UBOUND(CS%<f>), source=CS%<f>)`. Guard a plain
  `allocatable` field with `allocated(CS%<f>)`; guard an `ALLOCABLE_` field with the condition
  init allocates it under, never `allocated()` (the macro is empty under static memory, where
  `allocated()` doesn't compile).
- **Free** every member in the type's `*_end` (`%free()` is safe when unallocated).
- **Readers** take the bundle and `%view` each member they use into a local pointer named after
  the field, nullified first and viewed under `%associated()`, so `CS%<f>(i,j)` becomes
  `<f>(i,j)` and the arithmetic is unchanged. **Name the dummy differently from its type**
  (`metrics` for `hv_metrics_CS`): a variable may not share a name with an accessible derived
  type (gfortran: "Derived type 'hv_metrics_cs' … is being used before it is defined"). A
  component may share it.
- A field the tree writes is state, not configuration: leave it raw for
  `convert_array_containers`.
- The mirror and converter follow Step 4's container rule.

### 2. Trace usage across the call tree

For every subroutine from the entry point down, record which eligible fields it reads, writes or
forwards (via flattened dummies, `<dummy>%<field>`, or the whole struct). Beyond a handful of
subroutines, build the map with an Explore/fork agent.

### 3. Cluster by co-occurrence

Group fields read or forwarded together across most of the same subroutines. Fields without a
recurring companion stay as found.

### 4. Define each cluster's bundle

In the original type's module:

```fortran
!> <One-line purpose of this cluster of options>.
type, public :: <cluster>_CS
  <field's own declared type> :: field1 !< <doc, copied from the original struct's field doc>
  ...
end type <cluster>_CS

!> bind(C) mirror of <cluster>_CS, field-for-field, same order.
type, bind(C) :: <cluster>_CS_C
  <c_kind or _C container-mirror type> :: field1
  ...
end type <cluster>_CS_C

!> Converts a <cluster>_CS to its bind(C) mirror.
function <cluster>_CS_to_c(opts) result(cdesc)
  type(<cluster>_CS), intent(in) :: opts
  type(<cluster>_CS_C) :: cdesc
  cdesc%field1 = opts%field1          ! scalar: direct assignment
  cdesc%field2 = opts%field2_a%to_c() ! container: delegate
  ...
end function <cluster>_CS_to_c
```

The converter is unconditional if every member is a scalar (`reconstruction_CS_to_c`). With any
container member it is wrapped in `#ifdef _TIM` (`%to_c()` exists only under TIM), with its doc
comment before the `#ifdef` (`BT_cont_container_to_c`).

C++ side: a scalar-only mirror needs no marshalling helper; the bridge passes
`const <cluster>_CS_C*` and the kernel reads `const <cluster>_CS_C&`. It is defined in the
module's kernel header, not `turbotmp_bridge_c_types.h` (protoMOMxx `generate_amrex_code` lessons.md
§2.2). Field order and kinds (`c_int`, `c_bool`, `c_double`) must match that header exactly.

### 5. Restructure the original type

Replace each cluster's fields with one `type(<cluster>_CS) :: <cluster>_CS` member. Every
non-clustered field stays where it is.

### 6. Per leaf, decide once: keeps the whole struct, or stops taking it

Every subroutine touching the type's fields, in any module, is exactly one of:

- **Keeps the whole struct** (no bridge now or planned): no signature change; only
  `<dummy>%<field>` → `<dummy>%<cluster>_CS%<field>` in the body.
- **Stops taking it** (bridged, or about to be): its dummies become the bundles it uses plus any
  unbundled field as a bare dummy (`convert_array_containers` Step 2), in the same pass.

Never both for the same type (the bundle is inside the struct). Decide for every leaf before
editing.

### 7. Roll out one cluster at a time

For one cluster, at each "stops taking it" leaf using its fields: replace whatever carries them
(flattened dummies or struct access) with the bundle dummy. A bridged leaf changes its
`_fortran`/shim/`bind(C)` triple together. Apply the hot-loop exception only at its one leaf.
Where a leaf had renamed a field's dummy, restore the bundle field's real name.

## Verify, by script

- Argument count, order and kind match across every touched `_fortran`/shim/`bind(C)` triple
  (or plain signature) and every call site.
- Each converter: scalars assigned directly, containers via `%to_c()`, `#ifdef _TIM` iff it has
  a container member.
- No bare reference to a bundled field name remains in the tree outside `%<field>` access.
- Every leaf is exactly one of the Step 6 kinds.
- Comments: every existing `!<` survives verbatim; line limit (lessons §7); `#if`/`#endif` depth-tracked across the
  file (never negative, ends at 0).
- **Build:** don't compile, build or run tests unless the user explicitly allows it in this
  conversation (`~/.claude/CLAUDE.md`). Report it as unrun, for both infra layers.

## Invariants (hard rules)

- Only a type private to the tree; never convert a raw array here (Step 1b aside); never bundle
  an un-bindable field.
- Clusters from traced usage; no bundle for a field without a recurring companion; hot-loop
  exception at one leaf only.
- Each leaf keeps the whole struct or takes bundles, never both.
- Never reword, shorten or "fix" an existing comment; no comments about the refactor itself.
- Never claim an unrun build passed; never install a compiler.
- Exactly one `!!SKILLS: 0.4.0` marker per touched Fortran file, on its own line after the
  license header, before `module`; update an existing one in place. The version is shared by
  the whole skill family and bumped in lockstep.

## Commit gating

As `convert_array_containers` Step 8: `--enable_git_commit`/`--disable_git_commit` override;
otherwise `~/.claude/preferences.json` `git_commit_and_push` (`"auto"` commits; anything else,
missing or unreadable, skips and lists the files). Branch:
`claude_<lowercased_type-name>_config_bundle`. Never force-push.

## Output to the user

1. Clusters and their fields; rolled out vs remaining, if staged.
2. Fields left unbundled and why; Step 1b bundles built, if any.
3. Where the hot-loop exception was applied; naming inconsistencies fixed.
4. Excluded fields that gate computation, and the leaves they keep on the whole struct.
5. Verify results and build status; committed, or the files to commit.
