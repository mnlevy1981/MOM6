---
name: convert_locals_to_containers
version: "0.3.2"
description: Convert a subroutine's own internal scratch-local arrays (never its dummy arguments) to RealArray_t/IntArray_t/LogicalArray_t containers, wherever every callee that consumes a given local already accepts a container for that argument -- eliminating the alloc-source-copy-free round trip currently needed to hand a freshly-computed local off to an already-converted callee. Accepts a comma-separated list of subroutine names so a whole call tree (or however much of it has had convert_array_containers run on its interfaces already) can be swept in one invocation, each subroutine verified independently. Distinct from convert_array_containers (which only ever touches dummy arguments, never locals) and from hoist_container_marshalling (which only reorders/merges alloc, free and copy2F calls for containers that already exist -- it never changes a variable's fundamental type). Use this once a callee a subroutine calls has already been converted and that subroutine still computes the matching argument into a raw local purely to copy it in and immediately throw it away.
user-invocable: true
argument-hint: <work-directory> <function-name>[,<function-name>...] [--enable_git_commit] [--disable_git_commit]
---

# Convert a subroutine's scratch locals to containers, not just its interface

## Why this exists

`convert_array_containers` converts a subroutine's *dummy* arguments
and stops there — it never looks at locals. Once a callee is
converted, a caller that computes that callee's argument into a plain
local (rather than one of its own dummies) is stuck doing
`alloc(source=local)`/call/`free` at every call site, purely to hand
off a value just computed and never read again raw — the pattern seen
in `zonal_mass_flux`/`meridional_mass_flux`'s scratch locals at their
calls into `present_uhbt_or_set_BT_cont`/`present_vhbt_or_set_BT_cont`.

A third, distinct kind of work in this family: `convert_array_containers`
decides raw-vs-container for **dummy arguments** and rewrites the
public interface; `hoist_container_marshalling` reduces
`%alloc`/`%free`/`%copy2F`/`%copy2Array` **churn for containers that
already exist**, never changing a variable's type; this skill decides
raw-vs-container for a **local, non-dummy variable**, touching the
signature not at all.

Natural order down a call tree: `convert_array_containers` → (optional)
`convert_optional_args_to_containers` → **this skill** (stop rebuilding
containers from scratch for values already computed here) →
`hoist_container_marshalling`. This skill can create new containers for
that last step to group, but never changes the subroutine's own dummy
list — its effect is strictly downward, one level at a time, like the
other two.

Assumes `array_container_lessons` (container API §4, naming §7) has
already been invoked this session.

## Scope

**In scope:** local array variables (never a dummy) where at least one
callee they're passed to already takes a container for that exact
argument.

**Eligibility is per-local, not per-subroutine** — same granularity
`hoist_container_marshalling` uses for its own entanglement check, just
applied to a different question. A subroutine can have some locals
convert and others stay raw in the same pass; eligibility depends only
on where that particular local is used (see "The core technique" for
the case of a local feeding both a converted and a still-raw callee —
common, not a blocker).

**Out of scope:** converting a subroutine's own dummies
(`convert_array_containers`'s job, assumed already run); reordering/
merging existing `%alloc`/`%free` calls unrelated to this promotion
(`hoist_container_marshalling`, afterward); touching any callee's
signature or body; rewriting loop bodies or changing numerical
results — the `%view` pointer keeps the local's original name, so the
body needs zero edits beyond the declaration; silently reinterpreting
an `!$omp target`/`map(...)` directive naming a candidate local (see
Hard rules — whether a `%view`-derived pointer maps to a device the
same way a plain array does is compiler-specific and unverifiable on a
machine with no compiler).

## The core technique

For every local array in a target subroutine (skip dummies entirely):

1. **Never called anywhere → not in scope.** Pure internal scratch
   gains nothing from containerizing. Leave it.

2. **Called into at least one already-converted callee → eligible.**
   Check every call site the local feeds against that callee's current
   dummy at the matching position: if it's already a container, this
   site needs the alloc-copy-free round trip today — convert the local
   so it can pass the container directly. If the callee's dummy is
   still raw, nothing changes at that site (or it still passes a
   `%view` pointer, if the caller's own dummy feeding it was itself
   already converted) — and this does **not** block converting the
   local for its other, already-converted call sites. A local converts
   as soon as *any* of its call sites benefits; passing it as a
   container to one callee and a `%view` pointer to another in the same
   subroutine needs no special handling (ordinary guard-free `%view`,
   lessons §4.3) — a local is always mandatory, never `optional`, so
   none of the disassociated-pointer/`nullify` machinery (lessons §9
   #16) applies here.

3. **Referenced by an `!$omp target`/`map(...)`/similar directive by
   name → stop and ask, don't guess.** Retyping the local's storage
   from a plain array to a container's backing store changes what the
   directive maps, and this can't be verified by inspection on a
   machine with no compiler. Flag every such local, the same way
   `convert_array_containers` Step 2 stops on every `optional` array
   dummy. **Not a hypothetical for this skill's flagship use** — the
   `zonal_mass_flux`/`meridional_mass_flux` locals named above are all
   currently named in exactly this kind of
   `!$omp target enter data map(alloc: ...)` directive, so expect this
   question on the very first real run.

4. **Otherwise → convert.** Rename `<name>` to `<name>_a` (container,
   same rank), preserving its plain trailing comment verbatim if it has
   one (locals carry plain `!` comments, not dummies' `!<` Doxygen
   form — the check here is "the wording survives," not a marker
   count). Declare a pointer under the local's original name, `%view`
   it once; the rest of the body needs zero edits, same as converting
   a dummy.

   **Allocate with shape only, no `source=`.** Unlike a call-site
   container built from an externally-meaningful array (which defaults
   to copy-in, lessons §6, because the caller's raw array already held
   real data), a local being containerized never had a separate "old"
   value to lose. Confirm genuine write-before-read first (true of
   ordinary Fortran locals by default); a real read-before-write on
   some path needs copy-in instead, same as any `intent(inout)` case.

5. **Delete the old marshalling** at every call site that now receives
   the container directly — the local's own new container replaces the
   throwaway one that used to be built from it.

## Steps

### 0. Validate inputs

`$0` = work-directory, `$1` = comma-separated function-name list.

1. Help/empty-argument/flag handling identical to
   `convert_array_containers` Step 0 items 1–3.
2. Invoke `array_container_lessons` now if not already done this session.
3. Split `$1` on commas; locate each declaration (same grep as
   `convert_array_containers` Step 0 item 6); stop on 0 or >1 matches.
   A target need not be container-native itself — a fully raw
   subroutine can still have an in-scope local.
4. Drop any name with no local array declarations besides its dummies;
   report "nothing to do" for it without erroring the whole run.
5. **Name-collision check.** Before renaming `<name>` to `<name>_a`,
   confirm neither exists already in that subroutine. If one does, stop
   and ask — don't silently pick a different suffix.

### 1. Build the "before" table

Per target, per local: rank, every call site it feeds (callee,
position), whether that position is currently a container, and any
`!$omp`/similar directive naming it. Classify each as **eligible**,
**not called anywhere** (skip), or **flagged** (ask before proceeding).

### 2. Convert each eligible local

Per "The core technique" step 4: rename, add the `%view` pointer,
alloc shape-only (write-before-read verified first). Do this for every
eligible local in every target before touching any call site.

### 3. Update every call site the local feeds

Container callee → delete the old alloc/copy-in/free block, pass
`<name>_a` directly. Raw callee → leave the `%view` pointer `<name>`
passed unchanged.

### 4. Verify

- Rebuild the Step 1 table: every eligible local now has exactly one
  `%alloc`/`%free`; every container-callee site passes the container
  directly, not a rebuilt one.
- Plain-comment text preserved verbatim (not a `!<` count).
- Line length ≤ 100 on every added line.
- Diff review: no edits inside any loop body; no callee's signature or
  body touched; only local declarations and call-site marshalling changed.
- Every flagged (OMP-directive) local reported with the user's actual
  decision, not silently resolved.
- **Build:** don't compile, build or run tests yourself unless the user
  explicitly allows it in this conversation; by default that is the
  user's job (`~/.claude/CLAUDE.md`). Say plainly the build is unrun (both
  infra layers) and flag the OMP-directive interaction
  specifically as the top thing to confirm once possible — a genuinely
  new mechanism (a `%view` pointer standing in for a plain array under
  an existing device-mapping directive), not just the ordinary
  disassociated-pointer caveat.

## Versioning marker

Every Fortran file this skill creates or modifies gets a `!!SKILLS: 0.3.2`
marker line — the shared version for this whole skill family. If
missing, add it right after the file's license/header block, before
`module`; if present, update it in place. Grep-able
(`grep -rn "!!SKILLS:"`), meant to be stripped later.

## Hard rules

- Never skip or duplicate the `!!SKILLS: 0.3.2` marker.
- Never convert a local never passed to any call.
- Never let one still-raw callee block conversion of a local that also
  feeds an already-converted one — convert it, leave the raw site's
  `%view` pointer unchanged.
- Never silently convert a local named in an `!$omp target`/`map(...)`/
  similar directive — always stop and ask, no matter how many times in
  a row; it's not automatable without a compiler to verify the interaction.
- Never add a `source=` copy-in for a newly-containerized local without
  confirming a genuine read-before-write — default to shape-only, the
  opposite of call-site marshalling's default, and state why.
- Never touch a subroutine's own dummy declarations, or a callee's
  signature or body.
- Never rewrite loop math.
- Do not add comments about the conversion itself — same standing rule
  as the rest of this family.
- Do not attempt to install a compiler, or claim an unrun build passed.

## Commit gating

Same as the rest of the family: `--enable_git_commit`/
`--disable_git_commit` override; otherwise
`~/.claude/preferences.json`'s `git_commit_and_push` key. Branch:
`claude_<lowercased_first_function_name>_local_containers`, or
`..._and_N_more_local_containers` for a multi-name run.

## Output to the user on success

Report, per target subroutine, then a combined total:
1. Which locals converted, and which call sites now pass the container
   directly vs. still pass a `%view` pointer to a still-raw callee.
2. Which locals were left alone and why (never called, or flagged —
   with the user's decision for each flagged one).
3. Plain-comment preservation, line-length, diff-scope confirmation.
4. Build status, with the OMP-directive caveat called out if any local
   converted under an approved directive.
5. Whether committed, or modified files for manual commit.
6. A one-line pointer to `hoist_container_marshalling` as the natural
   next step if new containers were created.
