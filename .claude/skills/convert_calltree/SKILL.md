---
name: convert_calltree
version: "0.3.2"
description: Given a call-tree entry point (a subroutine whose external signature must stay fixed, like continuity()), survey every derived type and array/optional argument reachable beneath it, classify each against the existing container/bridge-readiness skills, get a human decision wherever a fixed rule can't resolve one, record the plan, then execute it end to end. Use this instead of re-deriving the survey-and-decide process by hand for each new entry point.
user-invocable: true
argument-hint: <work-directory> <entry-point-name> [--enable_git_commit] [--disable_git_commit]
---

# Convert a call-tree entry point to be bridge-ready, end to end

> **STANDING RULE — use the pre-written skills; never re-implement them (user, 2026-10-23).**
> - Run this skill by invoking it with the Skill tool, at the start of Phase 1, 2 and 3 and on
>   every resume. Reading `SKILL.md` with Read/cat is not invoking it.
> - Also invoke `array_container_lessons` (and `cpp_bridge_lessons` before Phase 3) with the
>   Skill tool before the first MOM6 edit of the session.
> - Every transformation that a sibling skill performs is done by invoking that skill with the
>   Skill tool and following its procedure. This covers containerizing, shadows, bundles,
>   optional/present conversion, locals, hoisting, box iterators and bridges. Never reproduce
>   one by hand, or by an ad-hoc script that edits the source.
> - The only code this skill authors itself is Stage 1 (Step 5a): rename + wrapper and the
>   restructuring sub-steps the plan records. Before starting such work, say in the hand-off which
>   plan item it is and that no sibling skill covers it. Do it by verbatim line copies proved by
>   script.
> - Before each stage, list the skills the plan assigns to it and confirm each was invoked.
>   Name every invocation in the stage hand-off. If a needed transformation has no skill, stop
>   and ask; don't improvise one.

## Why this exists

Making a call tree bridge-ready takes several different targeted
transformations — containerizing raw arrays, converting optional
dummies, shadowing or bundling a derived type, bridging each leaf — and
none of the skills that perform them decide *which* one applies to
*which* target, or in what order. That survey-and-decide work was done
by hand for `continuity()`; it doesn't scale to the ~10 more entry
points expected to need the same treatment.

**Why the entry point is the unit of work.** Its external signature is
the one thing that can't change — everything below it is free to be
restructured, because nothing below it is visible to the rest of MOM6.
This skill adds no new transformation; it adds the missing layer above
the existing ones: survey, classify, decide once with a human, record
the decisions, then execute them.

## The three phases

**Phase 1 (interactive)** — survey the tree, classify every target, ask
about the ones a fixed rule can't resolve, write a durable plan.
**Phase 2 (checkpointed)** — build every container, bundle, and shadow;
stabilize the tree. **Phase 3 (checkpointed)** — bridge every
subroutine, leaf to root, once the tree is stable. Both execution
phases share the same discipline: no decisions left to make, but a
commit/push/CI gate between every stage, not one unattended run.

Phase 1 alone is a complete, valid use of this skill (e.g. to review
before committing to the campaign); Phase 2 and Phase 3 are separate,
later invocations that read the plan back. Phase 3 assumes Phase 2 is
done — it bridges a tree that's already fully container/bundle-ready,
never a still-raw one.

## Scope

An **entry point** is a subroutine whose external signature must not
change; its descendants, transitively, are fair game. A **target** is
anything Phase 1 must classify: a derived type referenced anywhere in
the tree, a raw array dummy, an optional dummy.

**Wrapper / TreeRoot convention.** The entry point's own body needs
somewhere to attach entry-level logic (a shadow's build/copy-back, a
bundle's construction) without touching the tree below it. Step 1d
classifies which of three cases applies; Phase 2 (Step 5a) acts on it
before anything else:

- **Case 1 — a wrapper already exists**, any shape or name. Nothing to do.
- **Case 2 — the entry point is a bodiless alias for a distinctly-named
  implementation** (e.g. a `use ..., only : entry=>impl` rename). Keep
  the implementation's name; author a new subroutine under the entry
  point's own name whose body just calls the implementation, and fix up
  the `use`/`public` lines.
- **Case 3 — no separate name exists at all** (one subroutine, called
  directly). Rename it `<entry>_TR` ("TreeRoot"), then author the
  wrapper the same way as Case 2, calling `_TR`.

Worked example: pre-campaign `continuity()` was a zero-body
`use MOM_continuity_PPM, only : continuity=>continuity_PPM` alias over
`continuity_PPM` — Case 2. `_TR` is a distinct convention from
`generate_cpp_bridge`'s own `_fortran` rename: `_TR` marks the tree's
root, once; `_fortran` marks a bridged leaf's original implementation,
per leaf, later.

Whichever case, the wrapper's body starts as a pure pass-through and
never gets bridged — every later stage that needs entry-level logic
edits it, never the implementation/`_TR` subroutine.

## Hard precondition

1. Confirm the entry point's signature really is externally fixed — if
   it isn't, there may be no reason to scope the campaign this tightly.
2. Check every descendant for callers *outside* this tree, repo-wide. A
   descendant reached from more than one entry point can't have its
   containerization decided in isolation — flag it in the plan (Step 4)
   as "shared with `<other entry point>`, decide jointly."
3. Find every **sibling call tree** (Step 1f): every other entry point in
   the host module that computes with the same control structure or
   shadowed types. Whether a test configuration calls it does not matter.
   Siblings are converted in the same campaign by default, because they
   share the structures each stage restructures. Leaving one out costs a
   second full campaign later, so only the user can exclude one (Step 3).

## Phase 1 — survey, classify, decide, record

### Step 1. Inventory

1a. The entry point's own current public dummy list — the contract
    Phase 2 must never change.
1b. Every derived type referenced (directly, or via `%field`), every
    raw array dummy, and every optional dummy, down to the leaves. Use
    an Explore/fork agent for anything beyond a handful of subroutines
    — don't re-derive a survey this size from memory.
1c. Every descendant's caller list, repo-wide (the shared-descendant
    check, above) — do this for the whole tree now, not per-target later.
1d. Which Scope-section case applies (1/2/3). Record the tree's actual
    root name, and whether Phase 2 must author a wrapper (Cases 2/3)
    and/or a rename (Case 3 only).
1e. Every runtime option (`get_param` flag, scheme selector, block
    size) that selects a different code path inside the tree, and
    which of those paths the test coverage actually exercises. Test
    coverage means the standard test cases **and the whole CI**, not
    just `double_gyre`. MOM6's CI (`.github/workflows/verify-*.yml`,
    `.testing/`) runs its own `tc*` configurations under symmetric and
    asymmetric memory, layout, restart, rotation, dimension-scaling,
    REPRO-vs-DEBUG, OpenMP and NaN tests, plus the unit tests, on
    Linux and macOS. Read the workflow files and the `tc*` parameter
    files to see which paths they reach; don't assume. The paths
    nothing exercises become the plan's list of override runs for
    bit-for-bit checks (Step 5). Example: CorAdCalc's `double_gyre`
    only runs `SADOURNY75_ENERGY` + `KE_ARAKAWA`. The other 8
    `CORIOLIS_SCHEME` values, `KE_UP3` (with and without its limiter),
    `BOUND_CORIOLIS`, `CORIOLIS_EN_DIS`, `CORIOLIS_ADV_NKBLOCK` and
    symmetric memory all need explicit `MOM_override` runs.
    Coverage decides the testing plan only, never the campaign's scope
    (1f).
1f. **Sibling call trees.** List every procedure in the host module and
    sort each into one of four kinds:
    - **in the tree:** reachable from the entry point (compute the call
      graph by script, including type-bound calls);
    - **sibling call tree:** an entry point of its own (usually `public`,
      called from another module) that computes with the tree's control
      structure, its shadowed types or its grid fields;
    - **lifecycle:** `*_init`, `*_end`, and scalar parameter updates such
      as `updateCFLtruncationValue`. These stay Fortran; record each one
      and why;
    - **plumbing:** converters, `record_*` helpers, `build_from`/
      `copy_back`.

    Then list every caller of each entry point and sibling, repo-wide
    (e.g. all four dynamics drivers), and note which runtime flag gates
    each call. Record every sibling in the plan (Step 4) as part of the
    campaign, with its own 1a–1e inventory. Never leave one out because
    no test configuration calls it: that is a testing problem, solved
    with an override run (1e), not a reason to shrink the scope.
    Example: the `vertvisc_family` survey left out `vertFPmix` because
    `FPMIX` is off in `double_gyre`. It shares `vertvisc_CS` with the
    three converted trees, so Stage 6a had to edit it anyway, and its
    conversion now needs a second campaign.

### Step 2. Classify each target against fixed rules

- Raw array dummy → `convert_array_containers`.
- Optional array/scalar dummy → `convert_present_to_associated`, after
  containerizing if it's an array.
- A struct/pointer only ever forwarded opaquely (confirm by grep, never
  assume from its name) → leave alone; record the confirmation.
- **Classify a shared or deferred struct by what the tree does with it,
  not by whether the test configuration allocates it.** If any
  subroutine in the tree writes it or computes with its fields, it is a
  container target (shadow or field containers), even when it is unused
  in the standard cases. Example: `vertvisc_family` first classified
  `ADp` as inert because `double_gyre` allocates none of its arrays, but
  `vertvisc` writes six of them, and Phase 2 needed an extra stage for
  it.
- A private, scalar-or-container-field control structure whose fields
  recur together across signatures → `create_config_bundle_type`.
- A struct used elsewhere in the repo outside this tree →
  `create_shadow_container_type`, not a wholesale conversion. "Used"
  means its *fields* are accessed outside the tree. A type with
  `private` components that outside code only declares and forwards
  (e.g. `CoriolisAdv_CS`, held by the dynamics drivers) is private to
  the tree.
- **`G` (`ocean_grid_type`) → fixed classification, no Step 3 measurement.**
  - Grid arrays come from the persistent `grid_core` / `grid_OBC` singletons
    (`src/core/MOM_grid_containers.F90`). The cap fetches them with `grid_core(G)` /
    `grid_OBC(G)` and passes `Gcore`/`Gobc` down (`convert_array_containers` Step 2b).
  - Grid index ranges reach the tree root as whole-column boxes that the cap builds from `G`,
    e.g. CorAdCalc's `bxH0 = [isc:iec, jsc:jec, 1:nz]` and
    `bxQ0 = [IscB:IecB, JscB:JecB, 1:nz]`. Staggered boxes are combined from their i/j ranges,
    never grown across staggerings.
  - `G` then stays only in the cap and in Fortran-only helpers the cap calls (e.g. diagnostics
    posting).
  - A field the tree needs that is in neither structure goes to the user.
- A **per-point kernel**: a `pure`/`elemental` procedure called once
  per grid point from inside a loop, taking scalars or fixed-size
  stencil arrays (`q4(4)`, `h8(8)`), never grid-shaped arrays (e.g.
  CorAdCalc's 23 WENO/UP3 helpers, continuity's `flux_elem`) → **not a
  container target and not a bridge target.** In Phase 3 it becomes an
  `AMREX_GPU_DEVICE AMREX_FORCE_INLINE` device primitive in the owning
  kernel's `<module>_kernel.hpp` (protoMOMxx `generate_amrex_code`, pointwise
  device primitive tier), called from inside that kernel's
  `ParallelFor`. A `bind(C)` call per grid point would be prohibitively
  slow. Leave it out of the Phase 3 wave computation; the box-level
  subroutine that calls it is the bridge target.

### Step 3. Quantify and decide what Step 2 couldn't resolve

Typically a large multi-field struct (`G`/`GV`/`US`-shaped) too big to
shadow or bundle without measuring first, or anything whose shared-vs-
private status isn't obvious. Also any proposal to leave a sibling call
tree (Step 1f) out of the campaign: present what it shares with the
tree and the cost of converting it later as a separate campaign, and
let the user decide. For each: measure total fields touched
vs. total fields, and whether the touched ones recur together (a
bundling/shadowing candidate) or are each used at only one or two sites
(not worth it). Present the measurement and the real tradeoff via
`AskUserQuestion` — one topic at a time, never bundled into one
multi-question call. Never pick a default silently.

### Step 4. Record the plan

Write a durable manifest to
`<work-directory>/.claude/calltree-plans/<entry-point-name>.md` — not
the ephemeral Plan-mode file. For every target: which skill handles it,
its settings (field lists, "leave alone" plus its grep confirmation,
"shared with `<X>`" flags), and the execution order below. This file is
what makes the campaign resumable without re-running Step 3, and what
Phase 2 and Phase 3 read.

**Record the surveyed base first:** the repo, branch and commit the
survey was done against, at the top of the plan. Every line number
and structural claim in the plan is only valid for that commit. (The
first CorAdCalc plan was surveyed against an older pinned commit, and
its loop-structure decisions turned out to be wrong for
`dev/turbo-debug`, which had since been rewritten.) Also record Step
1e's override-run list, and Step 1f's sorted procedure list with every
sibling call tree.

**List every bridge-readiness item as a Phase 2 stage.** For each
Phase 3 target, check it against the Phase 2 exit criterion (Step 6b)
and give every gap a Phase 2 stage: report or I/O splits, diagnostics
moved to the cap, `bind(C)` mirrors and converters for every container
group, shadows of written structs, field containers for whole structs.
Never defer Fortran-only work to Phase 3. Phase 3 should be only
bridges and C++.

**Structural restructuring belongs in Stage 1.** If Phase 1 decides the
tree must be restructured before containerization, record each change
as an ordered Stage 1 sub-step, each its own bit-for-bit checkpoint.
Examples: splitting a dispatch body into per-scheme subroutines;
replacing an in-kernel block loop with box/tile iteration; making a
per-point helper family `pure` and its caller loops `do concurrent`.
CorAdCalc ran 1a+1b (rename + wrapper, then box conversion with a
per-tile body subroutine), 1c (per-scheme split) and 1d
(`pure`/`do concurrent`). Put the Case-3 rename + wrapper (Step 5a) in
the first sub-step.

**Phase 2 execution order.** Each position is forced by some sub-
skill's own precondition or direction rule, not chosen arbitrarily;
record any deviation the survey justifies rather than reordering
silently.

1. **TreeRoot split** (Scope) — skip for Case 1.
2. **`create_shadow_container_type`**, per target classified "shared
   outside the tree" — built/copied-back inside the wrapper.
3. **`create_config_bundle_type`**, per private clustering candidate —
   skip if none.
4. **Optional-array containerization**, per subroutine: default to
   `convert_array_containers`; use `convert_optional_args_to_containers`
   instead, scoped to one pivot, when converting its children first
   would branch combinatorially on `present()`. Never both, for the
   same argument. `convert_present_to_associated` (item 8) is separate
   and later — not an alternative here. **Never leave an optional array
   dummy raw on a target this plan also schedules for Phase 3
   bridging** — `generate_cpp_bridge` requires every array dummy to
   already be a container (its own Step 0 precondition), so "leave as
   raw" is only a valid choice here for a target that will never be
   bridged. Check Phase 3's list (below) before answering this question.
5. **`convert_array_containers` — downward pass**, root to leaves,
   mandatory dummies (that skill's preferred direction for dummy
   conversion).
6. **`convert_array_containers` — upward pass**, leaves to root,
   `G`/`GV`/`US`-drop decisions and Step 2b promotions (that skill's
   *required* direction for these — undecidable top-down).
7. **`convert_locals_to_containers`**, per subroutine once its dummies
   and its callees' are stable — scratch locals feeding an
   already-converted callee.
8. **`convert_present_to_associated`**, per dummy about to cross a
   `bind(C)` boundary, ahead of Phase 3 bridging that leaf. Precondition:
   already a container (item 4).
9. **`hoist_container_marshalling`**, once, at the tree's root, last.
10. **`convert_loops_to_box_iterators`**, once the tree is stable and
    before Phase 3. It turns every scalar-range `do concurrent` into a
    loop over a core or derived `Box_t`, so each maps one-to-one onto a
    `ParallelFor`. Stage it by core box, one bit-for-bit checkpoint each
    (CorAdCalc steps 9c–9e).

**Phase 3 execution order.** Not a fixed list — computed from Step 1b's
call graph, restricted to every subroutine in the tree except the
wrapper (which never bridges) and the per-point kernels (Step 2, which
become device primitives rather than bridges; a subroutine calling only
per-point kernels counts as a leaf): a subroutine is ready to bridge once
every in-tree callee it still calls is already bridged, so leaves go
first and the tree's root goes last. Record the computed wave order
(the ready-together groups) in the plan; Phase 3 (Step 7) reads it.

## Phase 2 — checkpointed execution

### Step 5. Execute one stage at a time, gated by commit/push/CI

One branch for the whole run, `claude_<lowercased-entry-point>_calltree`,
checked out once before Stage 1 (confirm the tree is clean first). Pass
`--disable_git_commit` to every invoked sibling skill, unconditionally
— this skill's Step 5 is the only thing that commits, so the run lands
on one branch, not one per sibling.

**Before Stage 1 (and when resuming): check the plan is still current.**
Compare the plan's recorded base commit with the branch you are about
to change. If they differ, spot-check the plan's line anchors and
structural claims (loop structure, dispatch shape, signatures) against
the current source. If any no longer hold, return to Phase 1 for the
affected sections before executing anything.

**Execute the plan's decisions; don't re-ask them.** For each target,
read its recorded decision first (skill, settings, exclusions,
privacy/sharing status). Ask the user only about what the plan
explicitly left to the sub-skill, and apply that sub-skill's own rules
before asking. For example, CorAdCalc's plan fixed `CoriolisAdv_CS`'s
privacy, field scope and exclusions, and left only the exact
clustering to `create_config_bundle_type`. That skill's usage trace
settles the clustering, so no further question was needed.

For each of the 9 stages (Step 4's Phase 2 order, including any
Stage 1 sub-steps as separate checkpoints), in sequence:

1. **Do the work** — author code directly (Stage 1) or invoke every
   sibling skill the stage's targets are recorded against, each with
   `--disable_git_commit`. A stage often means several invocations; the
   checkpoint is per stage, not per invocation. "Invoke" means the Skill
   tool, never a hand-written equivalent (see the standing rule at the
   top). The hand-off lists every skill invoked in the stage.
2. **Verify** — every invoked skill's own Verify section, plus this
   skill's Step 6 checks scoped to just this stage's files (the
   full-tree sweep is Stage 9's job). If verification finds a problem,
   fix it before handing off. Run scripted static checks. A by-eye read
   is not enough:
   - every changed subroutine's arguments are all declared, with no
     undeclared uses;
   - no case-insensitive duplicate names (Fortran ignores case: a new
     `ke` clashes with an existing `KE`);
   - every call's argument count and order match its subroutine;
   - every executable line of moved code is present exactly once.
3. **Build and bit-for-bit test — the user's job.** Don't compile,
   build or run tests yourself, not even a single-file syntax check,
   unless the user explicitly allows it in this conversation (user rule,
   `~/.claude/CLAUDE.md`). Otherwise hand the build and the runs to the
   user, and never call an unrun build passed. Every stage must be
   bit-for-bit identical to the previous stage's output: the default
   test cases, plus the Step 1e override runs covering the code paths
   this stage touched. Name those runs explicitly in the hand-off,
   since the default cases may exercise only one path. A stage also
   isn't verified until the whole CI passes (item 4). Its rotation,
   dimension-scaling, layout, restart and REPRO-vs-DEBUG tests catch
   errors a bit-for-bit `double_gyre` run cannot.
4. **Commit and push** onto the branch, message naming the stage and
   what ran, if commit gating (below) allows. Otherwise hand the
   uncommitted diff to the user, who commits and pushes. Then check CI
   (confirm how this repo's CI is wired; don't assume). CI does not
   cover GPU offload paths, so say so whenever a stage changes
   OpenMP-offload or `do concurrent` behaviour.
5. **Stop and report** — stage, commit (or "uncommitted, for the
   user"), build/bit-for-bit/CI status — then wait. Never start the next
   stage in the same turn, and never proceed past a failed stage or an
   unresolved verification problem; surface it and let the user decide.

If a stage surfaces a target Phase 1 didn't anticipate, stop and return
to Phase 1 for it rather than guessing.

**Step 5a — Stage 1's work**, per Step 1d's case: Case 1 → nothing,
commit skipped. Case 2 → author the wrapper under the entry point's
name, calling the implementation's existing name. Case 3 → rename the
implementation to `_TR` first, then author the wrapper. (See Scope for
what each case means.) Plus any structural-restructuring sub-steps the
plan recorded for Stage 1 (Step 4), each its own checkpoint. These are
the only pieces of code this skill authors itself rather than
dispatching to a sibling. When a restructuring moves code between
subroutines, copy each executable line verbatim and prove it by script.
Only dispatch lines (`if`/`elseif` chains replaced by a `select case`)
and index-translation lines may change. Record the as-built result in
the plan, including any deviation from what the plan anticipated and
why.

### Step 6. Whole-tree verification (Stage 9)

One full sweep, before Stage 9's commit, alongside
`hoist_container_marshalling` itself: argument count/order across the
whole tree, zero stray references anywhere to a field that should have
moved, doc-comment/line-length/`#ifdef` balance for every touched file,
and the entry point's external signature (1a) still byte-identical.
Also confirm Step 1d's case played out correctly (Case 2: implementation
kept its name; Case 3: `_TR` exists, pre-split body unchanged), and that
the wrapper's own signature matches 1a and carries no `bind(C)`
interface.

### Step 6b. Phase 2 exit: the bridge-readiness audit

**Phase 2's purpose is to make the Fortran bridge-ready, and Phase 2 is
not done until it is.** Finishing the planned stages does not count.
Bridge-ready means, for every Phase 3 target in the tree and in every
sibling call tree (Step 1f):

1. Every dummy can cross `bind(C)`:
   - a container (`RealArray_t`/`IntArray_t`/`LogicalArray_t`) or a
     `Box_t`;
   - a non-pointer `real`/`integer`/`logical` scalar;
   - `grid_core_type`;
   - a derived type with a `bind(C)` mirror and a converter in the
     module (container groups need these too, not only scalar bundles);
   - or `OBC` as a pointer, which crosses as a `c_ptr` that the C++
     aborts on when non-null.

   No whole `G`/`GV`/`US`, control structure, `tv`, `VarMix`, `Waves`
   or diagnostics struct, and no `character` or `optional` dummies.
2. No Fortran-only work in the body: no `post_*`,
   `query_averaging_enabled`, halo passes, checksums, point-accel
   reporting, `MOM_error`, file I/O, `G%`/`GV%`/`US%` references or
   `SZ*_` macros. No container `%alloc`/`%free` either: C++ cannot
   allocate a Fortran container, so a lazy allocation inside a target
   (e.g. `vertvisc_coef_TR`'s `a1_shelf_u/v`) moves to the wrapper.
3. Every `do concurrent` iterates over a `Box_t`.
4. The per-point kernels (Step 2) take no `tv`/`VarMix`/`G`/`GV`/`US`.

Run this as a scripted audit (`audit_bridge_ready.py` in this skill's
directory) after every stage from Stage 10 onward, and include its output
in each hand-off. Any gap becomes another Phase 2 stage, recorded in the
plan. Declare Phase 2 done only when the audit prints `BRIDGE-READY` and
the user has committed that stage.

## Phase 3 — leaf-to-root bridging

### Step 7. Bridge one wave at a time, leaves first, gated by commit/push/CI

Same branch as Phase 2 (`claude_<lowercased-entry-point>_calltree`,
already checked out); do not create a new one. Pass
`--disable_git_commit` to every `generate_cpp_bridge` invocation,
unconditionally — it has its own commit step, and this skill's Step 7
is what commits instead, same reasoning as Phase 2's Step 5.

**Precondition check, per target, before invoking `generate_cpp_bridge`
on it — do not assume Phase 2 left every target ready.** Confirm every
array dummy is already `type(RealArray_t)`/`type(IntArray_t)`, and the
iteration domain already a `type(Box_t)` (`generate_cpp_bridge`'s own
Step 0 item 5). If a target fails this — most likely a Step 4-item-4
decision that left an optional array dummy deliberately raw before this
target was scheduled for bridging — stop for that target: either
containerize it now (`convert_array_containers`, a Phase-2-shaped fix
applied here because Phase 2 already ran) or exclude the target from
this campaign's bridging and record why. Never force a container
conversion through unreviewed, and never skip the target silently.

For each wave (Step 4's Phase 3 order — every target whose in-tree
callees are already bridged), in sequence:

1. **Do the work** — invoke `generate_cpp_bridge` on every target in the
   wave, each with `--disable_git_commit`, after its precondition check
   passes.
2. **Verify** — each invocation's own Step 9 (the three-mode matrix,
   `cpp_bridge_lessons` §17); CAPTURE mode is the real bar here, since
   no AMReX C++ implementation exists yet (that's a separate, later
   skill's job — see Step 8).
3. **Commit** onto the branch, message naming the wave and which
   subroutines it bridged, if commit gating allows. Otherwise hand the
   diff to the user (as in Step 5).
4. **Push**, then check CI.
5. **Stop and report** — wave, commit, CI status — then wait. Same
   rules as Phase 2's Step 5: never start the next wave in the same
   turn, never proceed past a failed wave or a precondition-check
   failure.

The tree's root (Step 1d's recorded root — `continuity_PPM`-equivalent)
is the last wave; the wrapper is never bridged, in any wave.

### Step 8. Whole-tree bridging verification

After the last wave: confirm every non-wrapper subroutine in the tree
has a `_fortran`/shim/`bind(C)` triple, the wrapper still has none, and
every shim's public signature is unchanged from what Phase 2 left it
(bridging must not alter a signature — only rename and wrap). Report
the AMReX side (the C++ implementation behind each `bind(C)` interface)
as the explicit next deliverable this skill does not produce. That is
protoMOMxx's `generate_amrex_code` skill, one invocation per bridged
subroutine, with `generate_amrex_unit_test` for its tests. List the
per-point kernels (Step 2) each bridged kernel calls, since
`generate_amrex_code` ports them as device primitives in the same
invocation rather than as bridges of their own. Also check the protoMOMxx
checkout is recent enough for the bridge API the skills describe: for
example, `turbotmp::make_array4` taking the array's own `lb`, and
`IntA4Box` for `LogicalArray_C`.

## Versioning marker

Every Fortran file this skill creates or modifies — directly (Step 5a)
or via an invoked sibling, in Phase 2 or Phase 3 — gets a
`!!SKILLS: 0.3.2` marker line, the shared version for this whole skill
family. Deliberately grep-able (`grep -rn "!!SKILLS:"`) and meant to be
stripped later.

## Hard rules

- Never perform a sibling skill's transformation without invoking that
  skill with the Skill tool, and never run this skill from a read of its
  file instead of invoking it (standing rule at the top).
- Never skip the `!!SKILLS: 0.3.2` marker on a file Step 5a edits directly.
- Never skip the shared-descendant check (1c) or the sibling survey (1f).
- Never leave a sibling call tree out of the campaign because no test
  configuration calls it; only the user can exclude one (Step 3).
- Never let test coverage decide scope. Coverage decides the override
  runs (1e); the whole CI must pass at every stage.
- Never declare Phase 2 done before the bridge-readiness audit (Step 6b)
  prints `BRIDGE-READY`, and never defer a readiness gap to Phase 3.
- Never let Phase 2 or Phase 3 resolve a judgment call Phase 1 didn't
  record — return to Phase 1 instead.
- Never re-ask a decision Phase 1 already recorded; read the plan
  first (Step 5).
- Never execute against a plan whose recorded base commit no longer
  matches the source it describes, without first re-checking it (Step 5).
- Never schedule a per-point kernel (Step 2) as a bridge target.
- Never bypass an invoked sibling's own hard rules or preconditions —
  this skill only decides *which* skill and *when*, never *how*.
- Never default an ambiguous classification (Step 3) silently.
- Never deviate from Step 1d's case: no cosmetic rename of an existing
  (Case 1) wrapper, no renaming a Case 2 implementation, no skipping
  Step 5a for Case 2/3.
- Never bridge the wrapper, and never confuse its `_TR` rename with
  `generate_cpp_bridge`'s `_fortran` rename — different marks, different
  points in the pipeline.
- Never let an invoked sibling commit or branch during Phase 2 or
  Phase 3 — always pass `--disable_git_commit`.
- Never advance to the next stage or wave in the same turn, past one
  with an unresolved verification problem, or past a red CI run.
- Never invoke `generate_cpp_bridge` on a target whose precondition
  check (Step 7) hasn't passed — fix or exclude the target first, never
  force a container conversion through unreviewed.
- Never bridge a target before every in-tree callee it still calls is
  already bridged (Step 7's wave order) — leaves first, root last.
- Don't compile, build or run tests unless the user explicitly allows
  it in this conversation; by default that is the user's job. Do not
  claim an unrun build passed; never attempt to install a compiler.

## Commit gating

Phase 2 and Phase 3 each commit on one branch — the same branch,
created once before Phase 2's Stage 1, never per-sibling branches.
Whether either phase commits/pushes at all follows the usual chain:
`--enable_git_commit`/`--disable_git_commit` on `convert_calltree`
itself, else `~/.claude/preferences.json`'s `git_commit_and_push` key.
Once committing, every invoked sibling still gets `--disable_git_commit`
regardless (Hard rules). Applies to Phase 2/3 only — Phase 1 produces a
plan file, not code changes. A user's instruction during the session
(e.g. "don't commit, that's my job") overrides both and means manual
for the rest of the campaign: leave every stage's changes uncommitted,
report `git diff --stat`, and let the user build, test, commit and
push before the next stage starts.

## Output to the user on success

**Phase 1:** the plan file path, and a summary table (target → skill →
decision) for review before Phase 2 runs, plus Step 1f's table of every
procedure in the host module (in the tree / sibling / lifecycle /
plumbing).
**Phase 2, per stage:** the stage number, what ran, the commit hash (or
"uncommitted" plus `git diff --stat`, under manual commit gating), the
build/bit-for-bit status with the override runs it needs, and CI status.
Then stop.
**Phase 2, after Stage 9:** Step 6's verification results, confirmation
the entry point's signature is unchanged, and a rollup of all 9 stages'
commits.
**Phase 2, done:** the Step 6b audit output (`BRIDGE-READY`) for every
target in the tree and its siblings.
**Phase 3, per wave:** the wave number, which subroutines it bridged,
the commit hash, CI status — then stop.
**Phase 3, after the last wave:** Step 8's verification results, and the
explicit list of `bind(C)` interfaces still needing an AMReX-side
implementation.
