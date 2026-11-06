---
name: convert_calltree
version: "0.4.0"
description: Make a whole MOM6 call tree bridge-ready, starting from an entry point whose external signature must stay fixed (e.g. continuity()). Surveys every derived type, array and optional argument beneath it and its sibling call trees, assigns each target to the right container/bundle/shadow/bridge skill, asks the user wherever no fixed rule decides, records a durable plan, then executes it stage by stage (Phase 2) and bridges leaf to root (Phase 3). Use instead of re-deriving the survey-and-decide process by hand for each new entry point.
user-invocable: true
argument-hint: <work-directory> <entry-point-name> [--enable_git_commit] [--disable_git_commit]
---

# Convert a call-tree entry point to be bridge-ready, end to end

> **STANDING RULE — use the pre-written skills; never re-implement them.**
> - Invoke this skill with the Skill tool at the start of Phase 1, 2 and 3 and on every resume;
>   reading `SKILL.md` is not invoking it. Also invoke `array_container_lessons` (and
>   `cpp_bridge_lessons` before Phase 3) before the first MOM6 edit of the session.
> - Every transformation a sibling skill performs (containers, shadows, bundles,
>   optional/present, locals, hoisting, box iterators, bridges) is done by invoking that skill and
>   following its procedure. A script may apply a skill's steps mechanically, never more; the
>   skill's own verification runs on its output.
> - If a skill lacks something the campaign needs, stop and propose a general addition to that
>   skill (or this one) for the user to review, then apply it by re-invoking the skill. Never
>   fill the gap with one-off code.
> - The only code this skill authors itself is Stage 1 (Step 5a). Before such work, say in the
>   hand-off which plan item it is and that no sibling skill covers it.
> - Before each stage, list the skills the plan assigns to it; name every invocation in the
>   stage hand-off. A needed transformation with no skill is a question for the user.

This skill adds no transformation of its own. It adds the layer above the sibling skills:
survey, classify, decide once with a human, record, execute. The entry point is the unit of
work because its external signature is the one thing that can't change; everything below it is
invisible to the rest of MOM6 and free to be restructured.

## The three phases

- **Phase 1 (interactive):** survey the tree, classify every target, ask about what no fixed
  rule resolves, write a durable plan. A complete use of the skill on its own.
- **Phase 2 (checkpointed):** build every container, bundle and shadow until the tree passes
  the bridge-readiness audit (Step 6b).
- **Phase 3 (checkpointed):** bridge every subroutine, leaf to root, on the stable tree.

Phases 2 and 3 are separate, later invocations that read the plan back; they make no new
decisions, and every stage or wave is gated by build/test/commit/CI, never one unattended run.

## Scope

An **entry point** is a subroutine whose external signature must not change; its descendants
are fair game. A **target** is anything Phase 1 must classify: a derived type referenced in the
tree, a raw array dummy, an optional dummy.

**Wrapper / TreeRoot.** Entry-level logic (a shadow's build/copy-back, a bundle's construction,
building boxes from `G`) needs a home above the tree. Step 1d records which case applies; Stage 1
(Step 5a) acts on it:

| Case | Situation | Stage 1 action |
|---|---|---|
| 1 | a wrapper already exists, any name or shape | nothing; never rename it |
| 2 | the entry point is a bodiless alias for a distinctly named implementation (`use …, only : continuity=>continuity_PPM`) | keep the implementation's name; author a wrapper under the entry point's name that calls it; fix `use`/`public` |
| 3 | one subroutine, called directly | rename it `<entry>_TR` (TreeRoot), then author the wrapper as in Case 2 |

The wrapper starts as a pure pass-through, is never bridged, and is where every later stage puts
entry-level logic. `_TR` marks the tree's root once; it is unrelated to `generate_cpp_bridge`'s
per-leaf `_fortran` rename.

**Fortran-only work below a bridged subroutine (host context).** When work that must stay
Fortran (a halo exchange, a call into another module) sits inside a subroutine to be bridged and
can't move up to the wrapper without changing the order of computation:
- move it into a Fortran-only **host helper** (never bridged);
- the wrapper gathers the Fortran objects the helpers need (`G`, `GV`, `CS`, …) into a
  module-private context type of pointers and passes `c_loc(context)` down the tree as one
  `type(c_ptr), value, intent(in) :: host` dummy; the wrapper dummies it points at get `TARGET`,
  which changes no caller;
- each helper takes `host`, recovers the context with `c_f_pointer`, and points locals named
  like the original dummies at its members, so its body is unchanged;
- data the bridged caller already holds as containers or boxes is passed explicitly, not through
  the context;
- the C++ side treats `host` like `OBC`: it aborts when a helper would be needed, until that path
  has a C++ implementation or a call back into the helper.

Example: `horizontal_viscosity` has five `*_host` helpers (QG Leith, Leith+E and GME smoothing,
ZB2020 copy).

## Hard precondition

1. The entry point's signature really is externally fixed.
2. Every descendant's callers are checked repo-wide (1c). A descendant reached from another entry
   point is flagged "shared with `<other>`, decide jointly".
3. Every sibling call tree is found (1f) and converted in the same campaign by default; only the
   user can exclude one (Step 3).

## Phase 1 — survey, classify, decide, record

### Step 1. Inventory

1a. The entry point's current public dummy list: the contract Phase 2 never changes.

1b. Every derived type referenced (directly or via `%field`), raw array dummy and optional dummy,
    down to the leaves, plus the call graph (computed by script, including type-bound calls).
    Use an Explore/fork agent for anything beyond a handful of subroutines.

1c. Every descendant's caller list, repo-wide, for the whole tree at once.

1d. Which Scope case applies; the tree's actual root name; whether Stage 1 needs a wrapper
    (Cases 2/3) and a rename (Case 3).

1e. Every runtime option (`get_param` flag, scheme selector, block size) that selects a code path
    in the tree, and which paths are actually exercised: the standard test cases **and the whole
    CI** (read `.github/workflows/verify-*.yml` and the `.testing/` `tc*` parameter files — they
    cover symmetric/asymmetric memory, layout, restart, rotation, dimension scaling,
    REPRO-vs-DEBUG, OpenMP, NaN and unit tests; don't assume). Every unexercised path becomes an
    override run in the plan's bit-for-bit list. Coverage decides testing, never scope.

1f. **Sibling call trees.** Sort every procedure in the host module into:
    - **in the tree:** reachable from the entry point (1b's call graph);
    - **sibling call tree:** its own entry point (usually `public`, called from another module)
      that computes with the tree's control structure, shadowed types or grid fields;
    - **lifecycle:** `*_init`, `*_end`, scalar parameter updates; these stay Fortran (record why);
    - **plumbing:** converters, `record_*` helpers, `build_from`/`copy_back`.

    List every caller of each entry point and sibling repo-wide, with the runtime flag gating
    each call. Every sibling enters the plan with its own 1a–1e inventory, even if no test
    configuration calls it — that is solved by an override run, not by shrinking scope, since a
    sibling sharing the restructured structures gets edited anyway and otherwise needs a second
    campaign.

### Step 2. Classify each target against fixed rules

Classify by what the tree **does** with a target (confirmed by grep), never by its name or by
whether the test configuration allocates it.

| Target | Skill / decision |
|---|---|
| raw array dummy | `convert_array_containers` |
| optional array or scalar dummy | `convert_present_to_associated` (after containerizing, if an array) |
| struct/pointer only forwarded opaquely | leave alone; record the grep confirmation |
| struct the tree writes or computes with, even if unused in the standard cases | a container target (shadow or field containers) |
| private control structure whose scalar/container fields recur together across signatures | `create_config_bundle_type` |
| struct whose *fields* are accessed outside the tree (declaring and forwarding a type with `private` components doesn't count) | `create_shadow_container_type`, never wholesale conversion |
| `G` (`ocean_grid_type`) | fixed, no Step 3 measurement: arrays via `grid_core`/`grid_OBC`, index ranges via boxes, as in `convert_array_containers` Step 2. The wrapper builds whole-column boxes (`bxH0`, `bxQ0`) from `G`; staggered boxes are combined from their i/j ranges, never grown across staggerings. `G` remains only in the wrapper and the Fortran-only helpers it calls. A field in neither structure goes to the user. |
| per-point kernel: `pure`/`elemental`, called once per point inside a loop, scalar or fixed-size stencil arguments (`q4(4)`), never grid-shaped | neither a container nor a bridge target. In Phase 3 it becomes an `AMREX_GPU_DEVICE AMREX_FORCE_INLINE` device primitive in the owning kernel's `<module>_kernel.hpp` (protoMOMxx `generate_amrex_code`); a `bind(C)` call per point is far too slow. The box-level caller is the bridge target. |

### Step 3. Quantify and decide what Step 2 couldn't resolve

Typically a large multi-field struct (`G`/`GV`/`US`-shaped) or anything whose shared-vs-private
status isn't obvious; also any proposal to exclude a sibling call tree (present what it shares
and the cost of a later separate campaign). Measure fields touched vs. total, and whether the
touched ones recur together (bundle/shadow candidate) or appear at one or two sites (not worth
it). Present the measurement and the tradeoff via `AskUserQuestion`, one topic per call. Never
pick a default silently.

### Step 4. Record the plan

Write `<work-directory>/.claude/calltree-plans/<entry-point-name>.md` (not the ephemeral Plan-mode
file). It is what makes the campaign resumable and what Phases 2 and 3 read. It records:
- **the surveyed base** (repo, branch, commit) at the top: every line number and structural claim
  is valid only for that commit;
- per target: the skill, its settings (field lists, "leave alone" with grep confirmation,
  "shared with `<X>`" flags);
- every policy sibling skills will need, so they apply it instead of asking: the device-directive
  policy (e.g. `!$omp target` removed when containerized, since the bridge becomes the GPU path),
  each optional array's decision, each bridge target's confirmed signature;
- Step 1e's override runs and Step 1f's sorted procedure list;
- the Phase 2 and Phase 3 orders below.

**Every bridge-readiness gap is a Phase 2 stage.** Check each Phase 3 target against Step 6b and
give every gap a stage (report or I/O splits, diagnostics moved to the wrapper, `bind(C)` mirrors
and converters for every container group, shadows of written structs, field containers). Phase 3
is only bridges and C++.

**Structural restructuring belongs in Stage 1**, as ordered sub-steps, each its own bit-for-bit
checkpoint, with the Case 2/3 wrapper (and rename) first. Examples: splitting a dispatch body into
per-scheme subroutines; replacing an in-kernel block loop with box/tile iteration; making a
per-point helper family `pure` and its caller loops `do concurrent`.

**Phase 2 execution order.** Each position follows from a sibling skill's precondition or
direction rule; record any deviation the survey justifies.

1. **TreeRoot split** (Scope) — skip for Case 1.
2. **`create_shadow_container_type`**, per struct shared outside the tree; built and copied back
   in the wrapper.
3. **`create_config_bundle_type`**, per private clustering candidate.
4. **Optional-array containerization**, per subroutine: `convert_array_containers` by default;
   `convert_optional_args_to_containers` instead, scoped to one pivot, when converting children
   first would branch combinatorially on `present()`. Never both for one argument. Leaving an
   optional array raw is valid only for a target Phase 3 will never bridge (`generate_cpp_bridge`
   requires containers).
5. **`convert_array_containers`, downward pass**, root to leaves: mandatory dummies.
6. **`convert_array_containers`, upward pass**, leaves to root: its Step 2 replacements and
   Step 4 drops (undecidable top-down).
7. **`convert_locals_to_containers`**, per subroutine once its dummies and its callees' are
   stable.
8. **`convert_present_to_associated`**, per dummy that will cross `bind(C)`; already a container
   (item 4).
9. **`hoist_container_marshalling`**, once, at the tree's root.
10. **`convert_loops_to_box_iterators`**, once the tree is stable; one bit-for-bit checkpoint
    per core box. Serial loop nests nobody has converted go first through
    `convert_serial_loops_to_do_concurrent`, as their own checkpoint.
11. **`consolidate_iteration_boxes`**, cleanup before the readiness audit and Phase 3: base boxes
    at the root, one commented grow per local box, purpose names; one checkpoint per new root box,
    then per group of kernels for the local boxes left (its "Staging" section).

**Phase 3 execution order** is computed from 1b's call graph over every subroutine except the
wrapper and the per-point kernels (a subroutine calling only per-point kernels is a leaf): a
subroutine is ready once every in-tree callee it still calls is bridged. Record the resulting
waves; the root is last.

## Phase 2 — checkpointed execution

### Step 5. Execute one stage at a time

One branch for the whole campaign (Phases 2 and 3), `claude_<lowercased-entry-point>_calltree`,
checked out once on a clean tree before Stage 1. Every invoked sibling gets `--disable_git_commit`,
always; only this skill commits.

**Before Stage 1 and on every resume:** if the plan's base commit differs from the branch, spot-
check its line anchors and structural claims (loop structure, dispatch shape, signatures)
against the current source; anything stale goes back to Phase 1 first.

**Execute recorded decisions; don't re-ask them.** Read each target's decision first. Ask only
about what the plan explicitly left to a sub-skill, and only after that sub-skill's own rules fail
to settle it. A target or judgment call Phase 1 didn't anticipate sends that item back to Phase 1.

For each stage (the Phase 2 order, with Stage 1 sub-steps as separate checkpoints):

1. **Do the work:** Step 5a for Stage 1, otherwise invoke every sibling skill recorded for the
   stage's targets (Skill tool, `--disable_git_commit`). The checkpoint is per stage, not per
   invocation.
2. **Verify:** each invoked skill's own verification, plus scripted static checks on this
   stage's files (a by-eye read is not enough): every argument declared, no undeclared uses; no
   case-insensitive name clashes (`ke` vs `KE`); call argument count and order match; every
   moved executable line present exactly once. Fix problems before handing off. Every stage, also
   run the Step 6b audit and include its output; its gap count must never rise.
3. **Build and test — the user's job.** Don't compile, build or run tests, not even a single-file
   syntax check, unless the user explicitly allows it in this conversation (`~/.claude/CLAUDE.md`);
   never call an unrun build passed, never install a compiler. Each stage must be bit-for-bit
   identical to the previous one on the default cases plus the Step 1e override runs for the
   paths it touched (name them in the hand-off), and the whole CI must pass.
4. **Commit and push** if commit gating allows, naming the stage and what ran; otherwise hand the
   uncommitted diff (`git diff --stat`) to the user. Check CI (confirm how it is wired). CI does
   not cover GPU offload: say so whenever a stage changes OpenMP-offload or `do concurrent`
   behaviour.
5. **Stop and report** (stage, skills invoked, commit or "uncommitted", build/bit-for-bit/CI
   status) and wait. Never start the next stage in the same turn or pass a failed stage, an
   unresolved verification problem or a red CI run.

**Step 5a — Stage 1's work.** The Scope case's action (Case 1: nothing, no commit; Case 2:
author the wrapper; Case 3: rename to `_TR`, then author the wrapper), then the plan's
restructuring sub-steps, each its own checkpoint. Moved code is copied verbatim line by line and
proved by script; only dispatch lines (an `if`/`elseif` chain becoming `select case`) and
index-translation lines may change. Record the as-built result in the plan, with any deviation
and why.

### Step 6. Whole-tree verification

One full sweep at the hoisting stage (item 9), before its commit: argument count and order across
the whole tree; no stray reference anywhere to a field that should have moved; doc comments, line
length and `#if` balance in every touched file; the entry point's signature (1a) byte-identical;
Step 1d's case played out (Case 2: implementation kept its name; Case 3: `_TR` exists with the
pre-split body); the wrapper matches 1a and has no `bind(C)` interface.

### Step 6b. Phase 2 exit: the bridge-readiness audit

Phase 2 is done when the Fortran is bridge-ready, not when the planned stages are. Bridge-ready
means, for every Phase 3 target in the tree and in every sibling call tree:

1. Every dummy can cross `bind(C)`: a container or `Box_t`; a non-pointer `real`/`integer`/
   `logical` scalar; `grid_core_type`/`grid_OBC_type`; a derived type with a `bind(C)` mirror, a
   converter and a `record_*` capture helper in the module, each covering every field of the type
   (container groups too, not only scalar bundles); `OBC` as a pointer (a `c_ptr` the C++
   aborts on when non-null); or the host-context `c_ptr` (Scope), handled the same way. No whole
   `G`/`GV`/`US`, control structure, `tv`, `VarMix`, `Waves` or diagnostics struct, and no
   `character` or `optional` dummy.
2. No Fortran-only work in the body: no `post_*`, `query_averaging_enabled`, halo passes,
   checksums, point-accel reporting, `MOM_error`, file I/O, `G%`/`GV%`/`US%` references or `SZ*_`
   macros. No `%alloc`/`%free`/`%rebound` of a dummy or a component of one — C++ cannot allocate a
   Fortran container, so lazy allocation of persistent or bridge-crossing data moves to the
   wrapper. Local scratch containers inside the Fortran implementation are fine (e.g. allocated
   once per call and `%rebound` per tile): the C++ side does its own tiling and scratch.
3. Every `do concurrent` iterates over a `Box_t`, except a loop whose range is clipped by OBC
   segment data: it runs only with OBC active, which the C++ side aborts on (item 1).
4. Per-point kernels take no `tv`/`VarMix`/`G`/`GV`/`US`.

Run it as the scripted audit `audit_bridge_ready.py` (this skill's directory). Every gap becomes
another Phase 2 stage in the plan; Phase 2 is done only when the audit prints `BRIDGE-READY` and
the user has committed that stage. Never defer a gap to Phase 3.

## Phase 3 — leaf-to-root bridging

### Step 7. Bridge one wave at a time, leaves first

Same branch as Phase 2. **Before invoking `generate_cpp_bridge` on a target**, confirm every array
dummy is a container and the iteration domain a `Box_t` (its Step 0). If not — typically an
optional array deliberately left raw — stop for that target: containerize it now with
`convert_array_containers`, or exclude it from bridging and record why. Never force a conversion
through unreviewed, never skip silently.

For each wave of the recorded order: invoke `generate_cpp_bridge` on every target in it
(`--disable_git_commit`); verify with each invocation's own three-mode matrix
(`cpp_bridge_lessons` §17), where CAPTURE mode is the real bar since no AMReX implementation
exists yet; then build/test, commit, CI, stop and report exactly as in Step 5 items 3–5. The root
is the last wave; the wrapper is never bridged.

### Step 8. Whole-tree bridging verification

After the last wave: every non-wrapper subroutine has its `_fortran`/shim/`bind(C)` triple, the
wrapper has none, and every shim's public signature is unchanged from Phase 2 (bridging only
renames and wraps). Report the AMReX side as the next deliverable, which this skill doesn't
produce: protoMOMxx's `generate_amrex_code` per bridged subroutine (with `generate_amrex_unit_test`),
listing the per-point kernels each one calls, since they are ported as device primitives in the
same invocation. Check the protoMOMxx checkout is recent enough for the bridge API the skills describe
(e.g. `turbotmp::make_array4` taking the array's own `lb`; `IntA4Box` for `LogicalArray_C`).

## Commit gating

Phase 1 writes only the plan. Phases 2 and 3 commit on the one campaign branch if
`--enable_git_commit`/`--disable_git_commit` on this skill, else `~/.claude/preferences.json`'s
`git_commit_and_push`, says so. A user instruction in the session ("don't commit, that's my job")
overrides both for the rest of the campaign: leave each stage uncommitted, report
`git diff --stat`, and let the user build, test, commit and push before the next stage.

## Invariants (hard rules)

- Sibling transformations only by invoking the sibling skill; never bypass its own preconditions
  or hard rules — this skill decides *which* and *when*, never *how*.
- The entry point's external signature never changes; the wrapper is never bridged; the Step 1d
  case is followed exactly (no cosmetic rename of a Case 1 wrapper or Case 2 implementation).
- Never skip the shared-descendant check (1c) or sibling survey (1f); only the user excludes a
  sibling; coverage never decides scope.
- No judgment call is defaulted silently, re-asked once recorded, or resolved outside Phase 1.
- Per-point kernels are never bridge targets; a target is bridged only after every in-tree callee
  it still calls.
- One stage or wave per turn, each bit-for-bit and CI-clean before the next.
- Never build, compile or run tests without explicit permission in this conversation.
- Exactly one `!!SKILLS: 0.4.0` marker in every Fortran file the campaign creates or modifies
  (directly or via a sibling), after the license header, updated in place; grep-able
  (`grep -rn "!!SKILLS:"`), to be stripped later.

## Output to the user

- **Phase 1:** the plan path; a target → skill → decision table; Step 1f's procedure table.
- **Phase 2, per stage:** stage, skills invoked, commit (or `git diff --stat`), build/bit-for-bit
  status with the override runs needed, CI status; after item 9, Step 6's results; at the end,
  the Step 6b `BRIDGE-READY` output for every target and sibling.
- **Phase 3, per wave:** wave, subroutines bridged, commit, CI status; after the last wave,
  Step 8's results and the `bind(C)` interfaces still needing an AMReX implementation.
