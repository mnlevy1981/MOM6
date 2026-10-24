---
name: convert_locals_to_containers
version: "0.4.0"
description: Convert a subroutine's own scratch-local arrays (never its dummy arguments) to RealArray_t/IntArray_t/LogicalArray_t containers wherever at least one callee consuming the local already takes a container for that argument -- eliminating the alloc-source-copy-free round trip used to hand a freshly computed local to a converted callee. Accepts a comma-separated list of subroutines so a call tree can be swept in one run, each verified independently. Distinct from convert_array_containers (dummies only) and hoist_container_marshalling (reorders marshalling of existing containers, never changes a variable's type).
user-invocable: true
argument-hint: <work-directory> <function-name>[,<function-name>...] [--enable_git_commit] [--disable_git_commit]
---

# Convert a subroutine's scratch locals to containers

`convert_array_containers` converts dummies only. A caller that computes a converted callee's
argument into a raw local must then `alloc(source=local)`/call/`free` at every call. This skill
makes the local itself a container. Order down a tree: `convert_array_containers` → (optionally)
`convert_optional_args_to_containers` → **this skill** → `hoist_container_marshalling`, which
groups the containers this one creates. Container API and naming: `array_container_lessons`
(§4, §7); invoke it once per session first.

**Scope.** Local arrays only, decided per local. Never a dummy, a callee's signature or body, or
loop math (the `%view` pointer keeps the local's name, so the body is untouched). Reordering
unrelated marshalling is `hoist_container_marshalling`'s job.

**Campaign decisions come from the caller.** Decisions recorded in the calling plan (e.g. a
device-directive policy) are applied, not re-asked; ask only where none exists.

## Help message

If the arguments are empty, `help`, `--help` or `-h`, print this and stop:

```
Usage: /convert_locals_to_containers <work-directory> <function-name>[,<function-name>...] [--enable_git_commit] [--disable_git_commit]

Convert scratch-local arrays (never dummies) of one or more MOM6 subroutines to
containers wherever a callee they feed already takes a container, removing the
alloc/copy/free round trip at those call sites.

Arguments:
  <work-directory>       An existing TURBO-ESM/MOM6 checkout (src/, config_src/).
  <function-name>,...    Comma-separated subroutines to process (case-insensitive).
  --enable_git_commit    Commit and push this run (overrides the global preference).
  --disable_git_commit   Don't commit this run (overrides the global preference).
```

## The core technique

For each local array:

1. **Passed to no call → leave it.** Internal scratch gains nothing.
2. **Passed to at least one container dummy → eligible.** One such site is enough; a site whose
   callee is still raw doesn't block it and simply gets the `%view` pointer. A local is never
   optional, so a plain `%view` suffices (lessons §4.3) and no `nullify` machinery (lessons §9
   #16) applies.
3. **Named in an `!$omp target`/`map(...)` or similar directive → flagged.** Retyping its storage
   changes what the directive maps, which inspection can't verify. Apply the plan's recorded
   device-directive policy; with none, stop and ask, per local.
4. **Otherwise → convert.** Rename to `<name>_a` (same rank), keeping its trailing `!` comment
   verbatim; declare a pointer under the original name and `%view` it once. Allocate with
   **shape only, no `source=`** — a local has no prior value — after confirming write-before-read;
   a genuine read-before-write on some path needs copy-in (lessons §6), stated in the report.
5. **Delete the old marshalling** at every site that now takes the container directly.

## Steps

### 0. Validate inputs

1. Arguments and commit flags as `convert_array_containers` Step 0 items 1–2.
2. `array_container_lessons` invoked this session; if not, invoke it now.
3. Split the name list on commas; locate each as `convert_array_containers` Step 0 item 4 (stop on
   0 or >1 matches). A target may itself be fully raw.
4. A target with no local arrays: report "nothing to do" for it and continue.
5. **Name collision:** if `<name>_a` already exists in the subroutine, stop and ask; never pick
   another suffix.

### 1. Build the "before" table

Per target, per local: rank, each call it feeds (callee, position, container or raw), and any
directive naming it. Classify as eligible, not called (skip), or flagged.

### 2. Convert each eligible local

Per core-technique item 4, for every eligible local in every target, before touching call sites.

### 3. Update every call site the local feeds

Container dummy → delete the alloc/copy-in/copy-back/free block, pass `<name>_a`. Raw dummy →
pass the `%view` pointer `<name>`, unchanged.

### 4. Verify, by script

- **Table:** every converted local has exactly one `%alloc` and one `%free`; every container site
  passes it directly.
- **Comments:** every local's plain-comment text survives verbatim; no comment line edited.
- **Diff scope:** only local declarations and call-site marshalling changed.
- **Format:** every added line within 100 columns.
- **Flagged locals:** each reported with the decision applied (plan or user).
- **Build:** don't compile, build or run tests unless the user explicitly allows it in this
  conversation (`~/.claude/CLAUDE.md`). Say the build is unrun, for both infra layers, and name a
  `%view` pointer under an existing device-mapping directive as the first thing to confirm.

## Invariants (hard rules)

- Only locals passed to at least one container dummy convert; a raw callee never blocks one.
- A directive-named local is never converted without a recorded policy or the user's answer.
- Shape-only allocation unless a read-before-write is shown.
- No dummy, callee, or loop math touched; no comment about the conversion added.
- Never claim an unrun build passed; never install a compiler.
- Exactly one `!!SKILLS: 0.4.0` marker per touched file, after the license header (update in
  place).

## Commit (gated)

As `convert_array_containers` Step 8, on branch
`claude_<lowercased first function-name>_local_containers` (`..._and_N_more_local_containers`
for several names); when skipped, list the modified files.

## Output to the user

Per target, then a total:
1. Locals converted; which sites now take the container and which still get a `%view` pointer.
2. Locals left alone and why, with the decision applied to each flagged one.
3. The Step 4 results and build status, with the device-directive caveat if it applies;
   committed, or the files to commit.
4. If containers were created: `hoist_container_marshalling` as the next step.
