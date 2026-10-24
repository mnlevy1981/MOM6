---
name: convert_loops_to_box_iterators
version: "0.4.0"
description: Convert the scalar-range `do concurrent` headers in a MOM6 subroutine (or a call tree) to iterate over `Box_t` boxes -- a few core boxes, one per staggering, built per tile by the tree root, plus derived boxes grown from them -- instead of scalar bounds such as `Isq:Ieq+1`. Every range is reduced to core-box-plus-offset and proven equal in both memory modes; loop bodies never change; unused scalar bounds and box dummies are dropped. Use after the tree has box/tile iteration and before Phase 3 bridging, so each loop maps onto one AMReX `ParallelFor(box, ...)`.
user-invocable: true
argument-hint: <work-directory> <function-name>[,<function-name>...] [--core <box>[,<box>...]] [--enable_git_commit] [--disable_git_commit]
---

# Convert scalar-range loops to box iterators

## Why this exists

A kernel that reads scalar bounds from its boxes (`Isq = bxQ%idxS(1)`, `array_container_lessons`
§5) keeps headers like `do concurrent (k=ksc:kec, j=Jsq:Jeq+1, I=Isq:Ieq)`. The AMReX port needs a
box per loop (`ParallelFor` over an `amrex::Box` built from a `Box_C`, protoMOMxx `generate_amrex_code`
lessons.md §3). The ranges in a tree reduce to a handful of **core boxes**, one per staggering
(`bxH` h/h, `bxU` B/h, `bxV` h/B, `bxQ` B/B, plus any scheme-specific range such as `bxQs`),
each grown by at most a cell or two. This skill makes that explicit: every loop iterates over a
core box or a derived box built once from one.

Campaign decisions recorded in the calling plan (`convert_calltree`) — which core boxes exist,
staging by core (`--core`) — are applied, not re-asked. Assumes `array_container_lessons` (§5
`Box_t` API and memory modes, §9 #3, §9 #18) is loaded this session.

## Help message

If the arguments are empty, `help`, `--help` or `-h`, print this and stop:

```
Usage: /convert_loops_to_box_iterators <work-directory> <function-name>[,<function-name>...] [--core <box>[,<box>...]] [--enable_git_commit] [--disable_git_commit]

Rewrite scalar-range do concurrent headers to iterate over core Box_t boxes
or boxes grown from them. Loop bodies are untouched.

  --core bxU,bxV         Convert only loops whose core box is one of these.
  --enable_git_commit / --disable_git_commit   Override the global commit preference.
```

## Scope

**In scope:** `do concurrent` headers whose every range is a box bound plus a constant offset.

**Out of scope (header left alone):** serial `do` loops (convert them first with
`convert_serial_loops_to_do_concurrent`); ranges clipped by non-box data (OBC
segment bounds, a `min`/`max`, a runtime stencil not already folded into a core box). Loop
bodies, index variable names/case/order, `DO_LOCALITY(...)` and `!$omp` directives never change.

## The core technique

1. **Tile-relative indices.** A loop over a tile-relative index (`kk=1:kmax`, because the scratch
   it uses is indexed that way) iterates a box in the same index space: the root sets the tile
   boxes' range for that dimension to `1:kmax` and passes the tile offset (`kstart`) to each
   kernel as a plain argument, so the body's `k = kstart + kk - 1` is unchanged. Never re-index
   scratch or rewrite bodies to fit a box.

   **Resolve every bound to its box** (via the extraction lines, the root's `%set` calls and the
   `G` fields behind them) and write each range as `core ± offset` per dimension.

2. **Prove it equal in both memory modes.** `IscB`/`JscB` are `isc`/`jsc` non-symmetric but
   `isc-1`/`jsc-1` symmetric, while `IecB = iec`, `JecB = jec` always. So an *end* offset may
   cross staggerings (`is-1:Ieq` is `bxH%growLo(1,1)`), a *start* offset may not (never grow
   `bxH` to reach `IscB` or shrink `bxQ` to reach `isc`): use the right core box, or add one.

3. **Add a missing core box in the tree root:** `safe_alloc` it beside the others, `%set` it per
   tile, free it after the tile loop (lessons §5). Box dummies lead every kernel's argument list
   in `bxH, bxQ, bxU, bxV, bxQs` relative order, with a `!<` comment giving the range; update
   every call site.

4. **Build each derived box once per subroutine**, never per loop: a local after
   `! Local variables` with a comment giving its range, built right after the extraction lines
   with one `grow`/`growLo`/`growHi`/`growBy` call (no temporaries), freed before
   `end subroutine` and before every early `return` (lessons §9 #3).
   ```fortran
   bxQ_pij = bxQ%growBy(lo=[0,0,0], hi=[1,1,0])
   ```
   Name it `<core>_<p|m|g><dims>` (`p` end, `m` start, `g` both), one group per differing
   dimension (`bxU_pi_mj`), a digit only for `n>1` (`bxQ_p2i`).

5. **Rewrite the header only**, keeping index variables, their case and order, and whether `k`
   is present (a 2-D loop uses dims 1–2 of a 3-D box):
   ```fortran
   do concurrent (k=bxQ_pij%idxS(3):bxQ_pij%idxE(3), j=bxQ_pij%idxS(2):bxQ_pij%idxE(2), &
                  I=bxQ_pij%idxS(1):bxQ_pij%idxE(1)) DO_LOCALITY(local(hArea_q))
   ```
   Past the line limit (lessons §7), move the `DO_LOCALITY` tail to a third line after `) &`. A loop that runs
   before the tile loop needs its core box allocated and set there (move the `safe_alloc` up,
   don't duplicate it); the tile loop re-`set`s it as before.

6. **Drop what became unused**, recomputing usage over the body: scalar bounds (extraction line
   and declaration) and box dummies (and the matching actual). Keep anything still read by an
   out-of-scope loop or an automatic array's bounds.

## Steps

### 0. Validate inputs

1. Arguments and flags as in `convert_array_containers` Step 0 items 1–2; `--core` limits the
   run to the named core boxes.
2. Invoke `array_container_lessons` if not yet loaded this session.
3. Each subroutine found exactly once; each already has box iteration (box dummies, or a tile
   loop that builds them), **or receives its bounds as integer dummies from a caller that has
   it**. Such a routine takes the caller's core boxes instead (technique 3, updating every call
   site) and its integer bound dummies are dropped once unused (technique 6). A bound the caller
   chooses at run time (e.g. `is_Kh` = `Isq` or `is-1` by option) is a scheme-specific core box
   the root sets from those scalars, never a derived box. Anything else: stop, it needs the
   box/tile conversion first.

### 1. Inventory and reduce

Table every `do concurrent` header: subroutine, line, raw range, reduced `core ± offset`; plus
the distinct ranges, core boxes needed, and out-of-scope loops with the reason. **Show it before
editing.** A range with no core box, or a new core box the plan doesn't record, is a question.

### 2. Name-clash check

Check every new name against every symbol in the targets and the module, ignoring case (lessons
§9 #18). On a clash, stop and ask.

### 3. Edit

Per subroutine, in technique order 3 → 4 → 5 → 6. A script is fine; match each header against
its exact Step 1 range text, never a looser regex.

### 4. Verify

Scripted:
- **Canonicalised diff:** with each header joined to one line and its ranges replaced by a
  placeholder, the only differences from the original are box declarations/builds/frees,
  dropped bounds, and argument lists/call sites. Zero body lines.
- **Ranges:** every converted range, re-derived from its box, equals the Step 1 raw range in
  **both** memory modes.
- **Completeness:** no in-scope header keeps a scalar range; every derived box freed on every
  exit; every dropped name gone everywhere (repo `CLAUDE.md` "Verification discipline").
- **Calls:** argument count and order match every callee.
- **Declarations** precede the first executable statement; added lines within the line limit (lessons §7).
- **Build:** don't compile, build or run tests unless the user explicitly allows it in this
  conversation. Say it is unrun; suggest double_gyre bit-for-bit under both infra layers plus
  runs reaching each converted loop (every scheme branch, `KE_SCHEME` variant, a stencil that
  changes `bxQs`).

## Versioning marker

Exactly one `!!SKILLS: 0.4.0` line per touched Fortran file, after the license header, before
`module`; update in place if present.

## Hard rules

- Loop bodies, index variables, their order and directives never change.
- No derived box whose *start* crosses staggerings; every range proven in both memory modes.
- Derived boxes built once per subroutine and freed on every exit path.
- No serial loop or non-box-clipped range converted.
- No unused bound or box dummy left, none dropped that is still read.
- Box dummies lead, in `bxH, bxQ, bxU, bxV, bxQs` order.
- No comments about the conversion; never claim an unrun build passed.

## Commit gating

As `convert_array_containers` Step 8. Branch `claude_<lowercased first name>_box_iterators`, or
`..._and_N_more_box_iterators` for several names.

## Output to the user on success

The Step 1 table with each derived box's build; out-of-scope loops; core boxes added, argument
lists changed, names dropped; Step 4 results and build status; committed, or the files to commit.
