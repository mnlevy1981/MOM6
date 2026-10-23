---
name: hoist_container_marshalling
version: "0.3.2"
description: Reduce redundant %alloc/%free/%copy2F/%copy2Array churn in a MOM6 subroutine that already marshals RealArray_t/IntArray_t containers around calls to converted descendants (a Case-A caller, in convert_array_containers terms), and group the surviving alloc calls before -- and free/copy-back calls after -- the block of descendant calls, so that block reads and profiles as pure computation rather than memory bookkeeping. Use this on a subroutine that already has containers. The eligibility check is per-container, not per-subroutine: a container whose call sequence touches a still-raw callee, or whose alloc/free is guarded by one of the subroutine's own optional raw array dummies, is entangled and stays exactly as it is -- but that only excludes that one container, not the whole subroutine, from hoisting; every other container is grouped and simplified normally (Step 1). An entangled container is named as a follow-up (convert_array_containers on the raw callee, or convert_optional_args_to_containers on the optional dummy), the same way convert_array_containers itself defers a hoisting candidate rather than solving it inline. It never converts a new raw array, and it never touches a callee's signature or body. Best applied once every container conversion feeding into this subroutine is settled and every call site inside it is visible in one pass -- exactly the moment convert_array_containers' own Step 8/10 "hoisting candidate" deferral points to.
user-invocable: true
argument-hint: <work-directory> <function-name> [--enable_git_commit] [--disable_git_commit]
---

# Hoist container marshalling out of a Case-A caller's call sequence

## Why this exists

`convert_array_containers` defers this work deliberately: its Step 8
emits a straightforward per-call-site `alloc`/`copy2F`/`free` block even
when a container's source is loop-invariant, and lists these as
**hoisting candidates** instead of fixing them inline — hoisting needs
to see every call site in the subroutine at once, not one callee at a
time.

The result, once several callees under one caller have each been
converted separately (the common case for a module entry point), is
container churn that was individually correct but collectively
wasteful: the same container freed and reallocated from an unchanging
source, or copied out and immediately copied into a different container
from that same raw array, purely to hand a value to the next call.

This skill is the deferred pass: given a subroutine whose container
conversions are done, find every such pattern and either eliminate it,
replace a free+realloc with an in-place `%copy2Array`, or recognize two
containers as one continuous value and merge them — then group what
survives so the subroutine's block of descendant calls reads as pure
computation, not memory bookkeeping.

This skill assumes `array_container_lessons` has already been invoked
this session (container API — in particular §4.4/§4.6 — and the
intent-mapping rules). Not reproduced here.

## Scope

**In scope:** the marshalling in one subroutine that already speaks
containers — reordering, merging, or eliminating `%alloc`/`%free`/
`%copy2F`/`%copy2Array` calls for every container not entangled (below),
plus the loop-invariant scalar computations that feed them (e.g.
hoisting `edge_h_min = 2.0 * GV%Angstrom_H` out of a branch that
computes it twice).

**Entanglement is checked per container, in Step 1 — not per
subroutine.** A still-raw callee or an optional raw dummy on this
subroutine's own signature disqualifies only the container(s) actually
touching it; both patterns are tagged, not refused:

- **A container's `%copy2F` feeds a callee that still takes a raw array
  at that position.** That specific `%copy2F` can never join the
  grouped end-of-sequence block (§4) — it stays immediately before the
  call it feeds. The same container's `%alloc` can still join the
  grouped top block if nothing else about it is entangled; only the
  pinned operation is excluded, not the container's whole lifecycle.
- **A container is built from one of this subroutine's own `optional`
  raw dummies**, so its `alloc`/`free` sits inside an
  `if (present(...))` guard (lessons §6). That pair can never move to
  the unconditional top/bottom blocks — doing so would make the
  container always "present" regardless of whether the caller actually
  supplied the argument (lessons §9 #8b).

Either pattern excludes only the entangled container (or, for the
first, just the pinned operation) from §1–§4 below — leave it exactly
where it is and name it in the Step 5 report as a follow-up
(`convert_array_containers` on the raw callee; `convert_optional_args_to_containers`
on the optional dummy — the latter often removes the entanglement
entirely rather than relocating it). Every other container in the same
subroutine hoists normally. Only if *every* container is entangled (or
already minimal) does this skill have nothing to do — report that and
name what would unblock each one, rather than refusing the whole
subroutine up front.

**Explicitly out of scope:**
- Converting any raw array dummy or local to a container
  (`convert_array_containers`'s job).
- Touching any callee's signature or body.
- Touching a `%view`-guarded optional-forwarding block (lessons §9 #16)
  beyond what this skill's own hoisting requires.
- Rewriting loop bodies, `do concurrent` bounds, or numerical logic — a
  hoist is numerically inert by construction (lessons §1); any
  difference afterward is a bug.
- Moving `cpu_clock_begin`/`cpu_clock_end` boundaries, or adding a new
  clock or other timing instrumentation, unless explicitly asked.
- Hoisting across an entanglement point identified in Step 1 — that's
  Step 1's classification to make, not a judgment call to revisit here.

## The core technique — one litmus test, applied three ways

Every eliminable pattern reduces to one question, asked at each
`%free()` followed (anywhere later in the subroutine) by a re-`%alloc()`
of the same container, and at each `%copy2F(X)` followed by a
*different* container's `%alloc(source=X)`:

> **Does anything else read or write the raw array `X` between these two
> operations?**

Walk the actual body to answer this — never assume from a variable's
name or position. If the answer is no, one of three fixes applies.
(Step 1 tags entangled operations first; this test then applies to
everything left untagged.)

### 1. Source never changes — hoist to a single alloc/free

If a container's source (`G%mask2dT`, a mandatory dummy this subroutine
never mutates, ...) is the same at every use site, allocate it once at
the top of the subroutine (or of whichever branch uses it, if the
branches are mutually exclusive at runtime — confirm this explicitly,
e.g. a `logical` set once near the top and never reassigned), pass the
same container to every call site that needs it, and free it once at
the bottom.

```fortran
! BEFORE -- freed and rebuilt from the same unchanging source, twice
call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
call zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, ...)
call mask2dT_a%free()
...
call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
call meridional_edge_thickness(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, ...)
call mask2dT_a%free()

! AFTER -- one alloc, one free, for the whole subroutine
call mask2dT_a%alloc(lb=LBOUND(G%mask2dT), ub=UBOUND(G%mask2dT), source=G%mask2dT)
...
call zonal_edge_thickness(bxC, h_in_a, h_W_a, h_E_a, mask2dT_a, ...)
...
call meridional_edge_thickness(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, ...)
...
call mask2dT_a%free()
```

A container written by a call (`intent(inout)`) hoists the same way, as
long as nothing *else* reads or writes its source between uses — see §2
and §3 for when it does need to change.

### 2. Source genuinely changes once — `%copy2Array` in place, not free+realloc

If a container's value must change partway through the subroutine, that
is not a reason to `%free()` and `%alloc()` again — `%copy2Array`
(lessons §4.4) overwrites an *already-allocated* container's payload
from a raw array of the same shape, with no deallocate/reallocate cycle:

```fortran
! BEFORE -- free, then realloc from a source that DID change once
call h_in_a%free()
call h_in_a%alloc(lb=LBOUND(h), ub=UBOUND(h), source=h)

! AFTER -- same container, updated in place
call h_in_a%copy2Array(h)
```

Requires the container's shape to be identical across every use —
confirm this from the declarations, not by assumption.

**Variant: alloc/free inside a loop where only the bounds move.** This is per-tile scratch
whose shape is fixed but whose index bounds track the tile, e.g. `lb(3) = k_start`.
- Allocate once before the loop for the largest iteration, then inside the loop call
  `%rebound(lb=…, ub=…)` and re-`%view` (lessons §4.5a). Free once after the loop.
- Any `!$omp target enter/exit data` for those arrays moves out of the loop with the
  alloc/free, **only if the user agrees**; otherwise leave the directives where they are.
- Worked example: `CorAdCalc_TR`'s 11 tile scratch containers.

### 3. Two containers are one continuous value — merge them

The deepest form of this pattern, and easiest to miss: if container `A`
is written by call 1 and then container `B` is built fresh from `A`'s
own copied-out value purely to feed call 2, ask whether call 2 could
simply take `A` directly. **Fortran matches actual arguments by type,
rank, and shape, not name.** If they match, `B` need not exist at all
for that use:

```fortran
! BEFORE -- A written, copied out to raw h, a DIFFERENT container B
! rebuilt from that same raw h, purely to feed the next call
call continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a, hin_a=hin_a)
call h_a%copy2F(h)
call h_in_a%copy2Array(h)
call meridional_edge_thickness(bxC, h_in_a, h_S_a, h_N_a, mask2dT_a, ...)

! AFTER -- h_a IS the value meridional_edge_thickness needs as its h_in;
! pass it directly, positionally, into that slot. No round trip at all.
call continuity_zonal_convergence(bxC, h_a, uh_a, dt, IareaT_a, hin_a=hin_a)
call meridional_edge_thickness(bxC, h_a, h_S_a, h_N_a, mask2dT_a, ...)
```

`B` may still be needed under its own name for a *different* use
elsewhere (sourced from a genuinely different raw array) — merging only
removes the uses where `A` and `B` would otherwise hold byte-identical
data. Do not delete `B`'s declaration or its other, still-needed
alloc/free if so.

This fix requires tracing the actual data lineage, not pattern-matching
on variable names — `h_a` and `h_in_a` looking similar is a coincidence
of this file's naming convention (lessons §7), not a signal.

## Steps

### 0. Validate inputs

`$0` = work-directory, `$1` = function-name.

1. If `$0` or `$1` is empty, or `$1` is `help`/`--help`/`-h` → print
   usage (mirroring `convert_array_containers --help`'s format) and stop.
2. If `$0` is not an existing directory → stop:
   `Error: work directory "$0" does not exist.`
3. Parse `--enable_git_commit`/`--disable_git_commit`; stop if both are
   passed, same wording as `convert_array_containers`.
4. If `array_container_lessons` has not been invoked this session,
   invoke it now.
5. Locate `$1`'s declaration under `$0/src` and `$0/config_src`. If it
   has no `type(RealArray_t)`/`type(IntArray_t)` dummies or locals at
   all, stop and say so — nothing to do on a still-raw subroutine.
6. Read the full subroutine and build the Step 1 table *before*
   changing anything. If every container is already at one
   `%alloc`/one `%free` with nothing tagged, stop and report there's
   nothing to hoist. If every container is instead tagged entangled,
   stop and name what would unblock each one — run that separately,
   then come back.

### 1. Build the "before" table

For every container variable declared in the subroutine, walk the body
in program order and record: each `%alloc` and its `source=`
expression, each `%free`, each `%copy2F`/`%copy2Array` and its
target/source expression, and which call(s) it feeds in between. A
short script (grep the extracted text for `%alloc(`, `%free()`,
`%copy2F(`, `%copy2Array(`, tally by container name) is the fastest way
to build this and to re-verify it afterward.

**Tag entanglement for every operation, in the same pass.** For every
`call <callee>(...)`, check whether the actual argument at each
position is a container or a raw array; a container's `%copy2F` feeding
a raw-array position gets tagged **pinned — still-raw callee**.
Separately, if this subroutine's own dummy list has an `optional` raw
array, tag the container(s) whose `alloc`/`free` sit inside its
`if (present(...))` guard **entangled — optional raw dummy**. Everything
else is untagged and eligible below.

This table makes every later step mechanical rather than judgment-based
— once container `X` is known to be alloc'd 4 times, all with
`source=Y`, and untagged, §1's fix applies without re-reading the
surrounding code.

### 2. Classify every branch that matters

If an `if`/`else` (or `select case`) runs exactly one of several
mutually-exclusive paths, fixed before any of them execute — confirm
this explicitly — treat the paths as symmetric for hoisting: a
container used in only one path, with an unchanging source, is still
hoistable to a single subroutine-wide alloc/free, since only one path
ever runs per call. Do not hoist across a branch boundary if the
condition could be re-evaluated or the paths aren't truly exclusive.

### 3. Apply the three fixes from the table

Skip anything tagged pinned/entangled in Step 1. For every remaining
untagged container, in order: unchanging source with multiple
alloc/free pairs → hoist (§1); source changes exactly once → `%copy2Array`
in place (§2); a second container built purely to receive the first's
output → merge (§3), tracing actual data lineage.

After every fix, a container's remaining `%alloc`/`%free` should each
appear exactly once, unless it is genuinely used in two disjoint,
mutually-exclusive scopes with different shapes, or you have a stated
reason it cannot be merged — state that reason in the Step 5 report
rather than leaving it unexplained.

### 4. Group the survivors

Move every surviving, untagged `%alloc` to a single block before the
subroutine's main branch (or before the descendant-call sequence, if
none). Move every surviving, untagged `%free`, plus any untagged
`%copy2F` with no reader between its current position and the end of
the subroutine, to a single block after that sequence — `%copy2F` calls
first, then `%free`.

Two legitimate reasons an operation stays mid-sequence: it's a `%copy2F`
feeding a later `%copy2Array`/`%alloc(source=...)` in the same sequence
(the §3 merge case, or a pairing a merge turned out not possible for);
or it's tagged pinned/entangled in Step 1. Either way, leave it exactly
where it is — do not move a copy-back earlier or later than the point
that makes it correct.

The goal: the subroutine's block of descendant calls reads as pure
computation for every untagged container, with no marshalling calls
interleaved. A pinned `%copy2F` feeding a still-raw callee is the one
remaining exception, called out explicitly in the Step 5 report.

### 5. Verify

- Rebuild the Step 1 table against the edited code. Every untagged
  container's `%alloc` count must equal its `%free` count.
- Every pinned/entangled operation is byte-identical to before, and
  named in the report with what would unblock it.
- For every remaining `%copy2F`/`%copy2Array`, state in one line why
  it's necessary — an unexplained one suggests a missed Step 3 fix.
- Diff review: no edits inside any loop body; no callee's signature or
  body touched; only the calling subroutine's marshalling and branch
  structure changed.
- Line length ≤ 100 on every added line (lessons §7).
- No doc comment (`!<`) added, removed, or reworded — this skill never
  touches a dummy declaration, so the count before/after must be
  identical, not just balanced.
- **Build:** don't compile, build or run tests yourself unless the user
  explicitly allows it in this conversation; by default that is the
  user's job (`~/.claude/CLAUDE.md`). Say plainly the build was not run
  and that it needs both infra layers. A hoist is numerically inert by construction, but that
  claim is only as good as the trace in Step 3.

## Versioning marker

Every Fortran file this skill creates or modifies gets a `!!SKILLS: 0.3.2`
marker line — the shared version number for this whole skill family,
not just this one skill. If the file doesn't already have one, add it
as its own line immediately after the license/header comment block,
before the `module` statement; if it already has one, update it rather
than adding a second line. Deliberately grep-able
(`grep -rn "!!SKILLS:"`) and meant to be stripped later.

## Hard rules

- Never skip the `!!SKILLS: 0.3.2` marker on a file this skill touches,
  and never add a second marker line if one already exists.
- Never hoist a container's `%alloc`/`%free`/`%copy2F` across an
  entanglement point tagged in Step 1 — leave it exactly where it is,
  and name what would unblock it, instead of refusing the whole
  subroutine over it.
- Never convert a raw array to a container here — that's
  `convert_array_containers`'s job, run separately.
- Never touch a callee's signature or body.
- Never move or add a `cpu_clock` boundary or other timing
  instrumentation unless explicitly asked.
- Never hoist an alloc/free across a point where something else — not
  itself — reads or writes its raw source, without confirming from the
  actual body that no such reader/writer exists. "Doesn't look like
  anything else uses it" is not a substitute for checking.
- Never merge two containers (§3) based on their names looking related
  — trace the value lineage through the body first.
- Never leave a `%copy2F`/`%copy2Array` call whose purpose can't be
  stated in one sentence; justify it in the Step 5 report or remove it.
- Do not add comments about the hoist itself anywhere in the code —
  this codebase documents ocean physics, not its own refactoring
  history (lessons §7). The Step 5 report is where that reasoning
  belongs, not the source.
- Do not attempt to install a Fortran compiler, and do not claim a build
  passed that was never run.

## Commit gating

Same as `convert_array_containers`: `--enable_git_commit`/
`--disable_git_commit` override; otherwise follow
`~/.claude/preferences.json`'s `git_commit_and_push` key (`"auto"` runs
the commit; `"manual"`, a missing file, a missing key, or unparseable
JSON all mean skip it and report modified files for manual commit).
Branch name: `claude_<lowercased_function-name>_hoist`.

## Output to the user on success

Report:
1. The subroutine, and the before/after `%alloc` count total (across all
   untagged containers) as the headline number.
2. Per untagged container: before → after count, and which fix (§1/§2/§3)
   applied, one line each.
3. Any untagged container left with more than one alloc/free pair, and why.
4. Every entangled/pinned container or operation found in Step 1, left
   unchanged, and what would unblock it.
5. Confirmation the descendant-call sequence is now uninterrupted by
   marshalling for every untagged container; if a `%copy2F` had to stay
   mid-sequence, say what it feeds.
6. Doc-comment count unchanged, line-length check, diff-scope confirmation.
7. Build status.
8. Whether committed, or the list of modified files for manual commit.
