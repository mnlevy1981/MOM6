---
name: hoist_container_marshalling
version: "0.4.0"
description: Reduce redundant %alloc/%free/%copy2F/%copy2Array churn in a MOM6 subroutine that already marshals RealArray_t/IntArray_t containers around calls to converted callees (a Case-A caller), then group the surviving allocs before, and copy-backs/frees after, the block of calls so it reads as pure computation. Eligibility is per container -- one tied to a still-raw callee or to an optional raw dummy of the subroutine stays where it is and is reported as a follow-up; the rest are hoisted. Never converts a raw array, never touches a callee. Use it on the "hoisting candidates" convert_array_containers Step 6 reports, once the conversions feeding the subroutine are settled.
user-invocable: true
argument-hint: <work-directory> <function-name> [--enable_git_commit] [--disable_git_commit]
---

# Hoist container marshalling out of a Case-A caller

`convert_array_containers` marshals per call site (`alloc`/call/`copy2F`/`free`) and reports
loop-invariant sources as hoisting candidates. Once several callees are converted, that leaves
the same container rebuilt from an unchanged source, or copied out only to be copied into another
container. This skill removes that churn, using the container API in `array_container_lessons`
(§4.4, §4.5a, §6); invoke that once per session first.

**Scope.** Only the marshalling of existing containers in one subroutine (plus loop-invariant
scalars feeding it, e.g. a duplicated `edge_h_min = 2.0*GV%Angstrom_H`). Not: converting any raw
array, touching a callee's signature or body, loop bodies or `do concurrent` bounds, a
`%view`-guarded optional-forwarding block (lessons §9 #16) beyond what the hoist needs, or
`cpu_clock` boundaries/instrumentation (unless asked). A hoist is numerically inert (lessons §1);
any difference afterward is a bug.

**Campaign decisions come from the caller.** Decisions recorded in the calling plan (e.g. a
device-directive policy) are applied, not re-asked; ask only where none exists.

## Help message

If the arguments are empty, `help`, `--help` or `-h`, print this and stop:

```
Usage: /hoist_container_marshalling <work-directory> <function-name> [--enable_git_commit] [--disable_git_commit]

Remove redundant %alloc/%free/%copy2F/%copy2Array churn in one MOM6 subroutine
that already marshals containers, and group the survivors around its block of
calls. Containers tied to a still-raw callee or an optional raw dummy are left
as they are and reported as follow-ups.

Arguments:
  <work-directory>       An existing TURBO-ESM/MOM6 checkout (src/, config_src/).
  <function-name>        The subroutine to clean up (case-insensitive).
  --enable_git_commit    Commit and push this run (overrides the global preference).
  --disable_git_commit   Don't commit this run (overrides the global preference).
```

## The core technique

At each `%free` later followed by a re-`%alloc` of the same container, and at each `%copy2F(X)`
followed by another container's `%alloc(source=X)` or `%copy2Array(X)`, ask:

> **Does anything else read or write the raw array `X` between the two operations?**

Answer it by walking the body, never from names or positions. If no, one of three fixes applies.

### 1. Source never changes — one alloc, one free

Allocate once at the top (of the subroutine, or of a branch if only that branch uses it — see
Step 2), pass the same container to every call, free once at the bottom. Holds for
`intent(inout)` containers too, as long as nothing else touches the source between uses.

### 2. Source changes — `%copy2Array` in place, not free+realloc

`call h_in_a%copy2Array(h)` overwrites an allocated container's payload (lessons §4.4). Requires
an identical shape at every use; confirm it from the declarations.

**Variant — alloc/free inside a loop where only the bounds move** (per-tile scratch, e.g.
`lb(3) = k_start`): allocate once before the loop for the largest iteration, `%rebound` and
re-`%view` inside it (lessons §4.5a), free once after. Matching `!$omp target enter/exit data`
moves out with them only under the caller's recorded device-directive policy or the user's
agreement; otherwise leave the directives where they are.

### 3. Two containers are one value — merge them

If container `A` is written by one call, copied out to raw `h`, and container `B` is rebuilt from
`h` only to feed the next call, pass `A` in `B`'s slot (actuals match by type, rank and shape,
not name):

```fortran
! before                                          ! after
call conv(bxC, h_a, uh_a, ...)                    call conv(bxC, h_a, uh_a, ...)
call h_a%copy2F(h)                                call meridional_edge_thickness(bxC, h_a, ...)
call h_in_a%copy2Array(h)
call meridional_edge_thickness(bxC, h_in_a, ...)
```

Trace the data lineage, never infer from similar names. Keep `B`'s declaration and its other
uses if it also holds a genuinely different value elsewhere.

## Steps

### 0. Validate inputs

Stop on the first failure with a one-line error.
1. Arguments and commit flags as `convert_array_containers` Step 0 items 1–2.
2. `array_container_lessons` invoked this session; if not, invoke it now.
3. Locate the declaration as `convert_array_containers` Step 0 item 4. No container dummy or
   local at all → stop, nothing to do.
4. Build the Step 1 table before editing. Every container already at one alloc/one free with
   nothing tagged → stop, nothing to hoist. Every container tagged → stop and name what would
   unblock each.

### 1. Build the "before" table, tagging entanglement

Per container, in program order: each `%alloc` (and `source=`), `%free`, `%copy2F`/`%copy2Array`
(and target/source), and the calls in between. A grep-and-tally script is the fastest way to
build it and re-check it later. Tag, per container rather than per subroutine:
- **pinned — still-raw callee:** a `%copy2F` feeding a raw-array position of a call. Only that
  operation stays put, immediately before the call; the container's alloc can still hoist.
  Follow-up: `convert_array_containers` on the callee.
- **entangled — optional raw dummy:** alloc/free inside an `if (present(...))` guard on one of
  this subroutine's own optional raw dummies. The pair never moves to the unconditional blocks
  (lessons §9 #8b). Follow-up: `convert_optional_args_to_containers` on that dummy.

Everything untagged is eligible.

### 2. Classify branches

Mutually-exclusive paths (`if`/`select case`) whose condition is fixed before any of them runs,
and never re-evaluated — confirm it — are symmetric: a container used in only one path, with an
unchanging source, can still take one subroutine-wide alloc/free. Otherwise don't hoist across
the branch.

### 3. Apply the fixes

For each untagged container: repeated alloc/free from an unchanged source → §1; changed source →
§2; a second container built only from the first's output → §3. Afterwards each has exactly one
`%alloc` and one `%free`, unless it lives in disjoint exclusive scopes with different shapes or
there is another stated reason, reported.

### 4. Group the survivors

All untagged `%alloc`s go in one block before the main branch (or the call sequence); all
untagged `%copy2F`s with no reader before the end, then all `%free`s, go in one block after it.
Only two things stay mid-sequence: a `%copy2F` feeding a later `%copy2Array`/`%alloc(source=)`,
and anything tagged in Step 1. Never move a copy-back away from the point that makes it correct.
Every early return must still reach the frees (lessons §9 #15).

### 5. Verify, by script

- **Table:** rebuilt on the edited code, alloc count = free count for every untagged container;
  every tagged operation byte-identical to before.
- **Copies:** every remaining `%copy2F`/`%copy2Array` has a one-line reason; one without suggests
  a missed Step 3 fix.
- **Diff scope:** no loop body, callee, or `cpu_clock` call touched; only marshalling and branch
  structure.
- **Comments:** `!<` count identical before/after; no comment line edited.
- **Format:** every added line within 100 columns (lessons §7).
- **Build:** don't compile, build or run tests unless the user explicitly allows it in this
  conversation (`~/.claude/CLAUDE.md`). Say the build is unrun, for both infra layers.

## Invariants (hard rules)

- Numerical result unchanged; no loop math edited.
- No hoist across a Step 1 tag, nor across an unverified reader/writer of the source; a tag
  excludes one container, never the whole subroutine.
- No merge (§3) without traced lineage.
- No raw array converted; no callee touched; no `cpu_clock` moved or added unless asked.
- No comment about the hoist added to the code.
- Never claim an unrun build passed; never install a compiler.
- Exactly one `!!SKILLS: 0.4.0` marker per touched file, after the license header (update in
  place).

## Commit (gated)

As `convert_array_containers` Step 8, on branch `claude_<lowercased function-name>_hoist`;
when skipped, list the modified files.

## Output to the user

1. Headline: total untagged `%alloc` count before → after; per container, before → after and the
   fix applied (§1/§2/§3); any left with more than one pair, and why.
2. Every tagged container or operation, left unchanged, and what would unblock it.
3. Whether the call block is free of marshalling; any `%copy2F` left mid-sequence and what it feeds.
4. The Step 5 results and build status; committed, or the files to commit.
