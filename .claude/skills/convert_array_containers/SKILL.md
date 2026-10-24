---
name: convert_array_containers
version: "0.4.0"
description: Convert a MOM6 Fortran subroutine from raw grid-shaped arrays to RealArray_t/IntArray_t/LogicalArray_t containers in place -- array dummies become containers, the body gets %view pointers so the math is untouched, and every call site is updated. Also finishes an already-converted subroutine: every remaining reference through a derived-type dummy (G%, GV%, US%, CS%, ...) becomes a persistent-container view, a box, or a plain scalar dummy, locals stop being sized by grid macros, and dummies the body no longer uses are dropped. Use when pushing the raw-array boundary up the MOM6 call tree. Does NOT rename anything, add a dispatcher shim, capture mode, or a bind(C) bridge -- those belong to generate_cpp_bridge.
user-invocable: true
argument-hint: <work-directory> <function-name> [--enable_src_validate] [--enable_git_commit] [--disable_git_commit]
---

# Convert a MOM6 subroutine to array containers

This skill is the **procedure**. Templates, the container API, intent rules, worked examples and
pitfalls live in `array_container_lessons` (§1–§10); invoke it once per session before this skill.

**Data structures only.** No rename, `*_fortran` variant, dispatcher, capture mode, `#ifdef _TIM`,
`%to_c` or `bind(C)` interface; those are `generate_cpp_bridge`'s job (lessons §10).

**Campaign decisions come from the caller.** When this skill runs inside a campaign
(`convert_calltree`), the plan's recorded decisions (signatures, optional-array choices,
device-directive policy, which groups exist) are applied as recorded, not re-asked. Anything the
plan does not cover and this skill does not decide is a question for the user; never improvise it.

**Scripts are allowed.** The steps may be applied by a script, provided the script only does what
the steps say and Step 7's checks run on its output.

## Help message

If the arguments are empty, `help`, `--help` or `-h`, print this and stop:

```
Usage: /convert_array_containers <work-directory> <function-name> [--enable_src_validate] [--enable_git_commit] [--disable_git_commit]

Convert one MOM6 subroutine from raw grid-shaped Fortran arrays to
RealArray_t / IntArray_t / LogicalArray_t containers, in place, or finish
one that is already converted (derived-type references, macro-sized locals,
unused dummies). The subroutine keeps its name; the body gets %view pointers
so the math is unchanged; every call site is updated.

Does NOT rename to *_fortran, add a dispatcher shim, capture mode, or a
bind(C) bridge -- use /generate_cpp_bridge for that, afterwards.

Arguments:
  <work-directory>       An existing TURBO-ESM/MOM6 checkout (src/, config_src/).
  <function-name>        The subroutine to convert (case-insensitive).
  --enable_src_validate  Also check the checkout is TURBO-ESM/MOM6 with a clean tree.
  --enable_git_commit    Commit and push this run (overrides the global preference).
  --disable_git_commit   Don't commit this run (overrides the global preference).

Without either commit flag, ~/.claude/preferences.json's "git_commit_and_push"
("auto" or "manual"; missing or unreadable means "manual") decides.

Example:
  /convert_array_containers /glade/derecho/scratch/sunjian/MOM6 zonal_mass_flux
```

## Step 0 — validate inputs

Stop on the first failure with a one-line error; never assume a default or create anything.
1. Both positional arguments present; the work directory exists (no cloning).
2. Flags recognized; `--enable_git_commit` and `--disable_git_commit` not both given.
3. `array_container_lessons` already invoked this session; if not, invoke it now.
4. Exactly one declaration of `<function-name>` under `src/` or `config_src/`; more than one, ask.
5. With `--enable_src_validate`: `git remote -v` names `TURBO-ESM/MOM6`, the tree is clean (never
   stash or discard), and report the current branch.
6. **Something to do.** Stop and report "already container-native" only if all of these hold:
   no raw array dummy, no reference through a derived-type dummy in the body (Step 2), no local
   sized by a grid macro (Step 3), and no unused dummy (Step 4). A `<name>_fortran` sibling means
   `<name>` is a bridge shim: convert it normally, but never edit inside `select case (mode)` or a
   `bind(C)` interface beyond renaming a reference; needing more is `generate_cpp_bridge`'s job.

## Procedure

### 1. Containerize the array dummies

For every dummy (lessons §6, §7):
- grid-shaped `real`/`integer`/`logical` array → `RealArray_t`/`IntArray_t`/`LogicalArray_t`,
  renamed `<name>_a`; written ones are `intent(inout)`, never `intent(out)` (lessons §9 #6);
- scalars, `logical` scalars, pointers (`OBC`) and derived types → unchanged here (Step 2 decides
  derived types).

An **`optional` array dummy** is converted or left raw by the caller's decision; if there is none,
ask, per dummy. A converted optional keeps `optional`, its `%view` is guarded by `present(<name>_a)`,
and every `present(<name>)` is renamed to match.

In the body, declare one pointer per container under the original name and rank
(`real, dimension(:,:,:), contiguous, pointer :: h`) and `call h_a%view(h)` before first use.
The body is otherwise unchanged. A pointer set only under a guard is `nullify`d before the first
guard (lessons §9 #16). A device directive (`!$omp target …`) in the routine follows the caller's
device-directive policy; with none recorded, ask.

### 2. Replace every reference through a derived-type dummy

Grep the body, including call arguments, for every `<dummy>%<field>` of every derived-type dummy
(`G`, `GV`, `US`, a control structure, …). Each distinct field becomes exactly one of:

| The field is … | It becomes … |
|---|---|
| an array with a persistent container copy (`grid_core`/`grid_OBC` for `G`; a bundle such as a metrics group for a control structure) | a `%view` of that copy: the subroutine takes the group (`Gcore`, `metrics`, …) and declares a pointer named after the field |
| an index range or extent (`G%isc`, `G%IecB`, `GV%ke`, …) | the matching bound of a `Box_t` dummy (`bxH0%idxS(1)`, `bxQ0%idxE(2)`, `bxH0%idxE(3)`, …), at any depth of the tree |
| any other scalar | a plain scalar dummy of the same name; callers pass `<their dummy>%<field>` |
| needed only by a Fortran-only callee (halo exchange, another module) | left for that callee: the subroutine is not a bridge target, or the call moves out (`convert_calltree` decides) |

Rules that apply to every row:
- Never a per-call `alloc(source=G%<f>)` or a per-field `<f>_a` dummy for data that has a
  persistent copy. An array with no persistent copy is a question for the user.
- If a local container already exists from the same source, promote it instead of adding a dummy.
- The only body edit is the literal substitution `<dummy>%<field>` → its replacement, in code,
  never in comments.
- Boxes are built from `G` once, by the caller that still holds `G` (normally the wrapper), and
  passed down as `intent(in)` dummies.

Decide this **bottom-up**: a subroutine's needs aren't known until its callees have settled
theirs (lessons §1). Containerizing dummies (Step 1) may go top-down.

### 3. Size locals from container dummies, not grid macros

A local (or a dummy left raw by decision) sized with `SZI_`/`SZJ_`/`SZK_`/`SZIB_`/`SZJB_` is
re-sized dimension by dimension from a **mandatory** container dummy with the same extent in that
dimension: `SZIB_(G)` → `u_a%lb(1):u_a%ub(1)`. Never from an `optional` dummy or a `%view` pointer.
A dimension with no such dummy stays as it is; report it. If the result exceeds the line limit,
change only line breaks (attribute form, `dimension(...) :: &`, keeps the comment intact).

### 4. Drop what the body no longer uses

After Steps 1–3, drop every dummy the body no longer references: `G`, `GV`, `US`, a control
structure, a box, anything. Keep one that is only forwarded to a callee that still needs it.

### 5. Rewrite the declarations and the header

- New dummies (`Gcore`, groups, boxes, scalars) get a new `!<` comment describing what they hold.
- Every existing `!<` comment survives word for word, including when a two-line declaration
  becomes one; wrap with `!!`, never shorten (lessons §7).
- Line limit (lessons §7): move `&` breaks, never drop or rename an argument to fit.

### 6. Update every call site

Find every call, including continuation forms. Build each call's new argument list from the new
dummy list; keyword arguments stay keyword arguments.
- **The caller is still raw (Case A):** marshal: allocate each container with `source=` (also for
  written arrays, unless the callee writes the whole extent), call, `copy2F` written ones, `free`
  (lessons §3, §6). An optional source needs a conditional allocation (lessons §9 #8b).
  Loop-invariant sources are reported as hoisting candidates, not restructured here.
- **The caller is already converted (Case B):** pass its own containers, groups, boxes and scalars
  straight through; delete any marshalling that has become redundant.
- Check every position: a container actual for a container dummy, a raw actual for a raw dummy
  (lessons §9 #9). Count and order alone don't catch a swap.
- Don't move `cpu_clock_begin`/`cpu_clock_end`; every early return must still reach the `free`s
  (lessons §9 #15).
- A still-raw callee gets the `%view` pointer, not the container; report it as a follow-up. An
  `elemental` callee keeps its signature (lessons §9 #10).

### 7. Verify, by script

- **Math unchanged:** every changed body statement differs from before only by a Step 1 rename
  or a Step 2 substitution.
- **References gone:** the Step 2 grep comes back empty; no grid macro sizes a local that has a
  mandatory container match.
- **Declarations:** no undeclared names, no case-insensitive duplicates or clashes with a type
  name, every dummy declared.
- **Calls:** count, order and container/raw kind match at every site.
- **Comments:** every removed `!<` text reappears verbatim; no comment line edited.
- **Guards:** every guarded `%view` pointer is nullified first; no `%associated()` inside a device
  construct (lessons §9 #17).
- **Format:** line limit; `#if` depth balanced.
- **Build:** don't compile, build or run tests unless the user explicitly allows it in this
  conversation (`~/.claude/CLAUDE.md`). Say the build is unrun, for both infra layers.

### 8. Commit (gated)

`--disable_git_commit` skips; `--enable_git_commit` runs; otherwise the global preference decides
(unreadable means skip). Commit on `claude_<lowercased function-name>_containers` with a message
naming the subroutine, its new containers and the call sites updated, ending with the attribution
trailer. If the push is rejected, stop and report; never force-push.

## Invariants (hard rules)

- The numerical result never changes; loop bodies change only by Step 1/Step 2 substitutions.
- Doc comments are preserved verbatim; no comment about the conversion itself is added.
- No subroutine rename and no bridge machinery.
- Written containers are `intent(inout)`.
- Persistent data is viewed, never copied per call.
- Optional arrays and anything else the plan doesn't decide are asked about, never assumed.
- Every call site is updated in the same change.
- Never claim an unrun build passed; never install a compiler.
- Exactly one `!!SKILLS: 0.4.0` marker per touched file, after the license header (update in place).

## Output to the user

1. The subroutine's new argument list, and what each change came from (Step 1–4).
2. Every call site updated, and any still-raw callee or hoisting candidate left as a follow-up.
3. Optional-array decisions; locals re-sized and any left as they were.
4. The Step 7 results and build status; committed, or the files to commit.
