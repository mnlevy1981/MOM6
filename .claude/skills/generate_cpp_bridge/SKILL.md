---
name: generate_cpp_bridge
version: "0.4.0"
description: Wrap an already-container-based MOM6 Fortran subroutine in a runtime-dispatched shim that selects, per kernel via an env var, the original Fortran code, a capture mode that records inputs and outputs for offline validation, or a C++/AMReX bridge called through bind(C). Use when porting a MOM6 kernel to AMReX while keeping the Fortran truth as the reference. Requires RealArray_t/IntArray_t dummies and a Box_t iteration domain -- run convert_array_containers first.
user-invocable: true
argument-hint: <work-directory> <function-name> [--enable_src_validate] [--enable_git_commit] [--disable_git_commit]
---

# Generate a C++ bridge for a MOM6 subroutine

This skill is the **procedure**. Templates, type maps, naming and pitfalls live in
`cpp_bridge_lessons` (§1–§17, cited here as "lessons §N"); invoke it once per session first.

**Bridge machinery only.** Converting data structures (array dummies → containers, `%view`
pointers, `Box_t` dummies and loops, call-site `alloc`/`copy2F`/`free`) is
`convert_array_containers`'s job. If the target needs any of it, stop and send the user there;
never do it here. This skill adds only the `*_fortran` rename, the shim, capture mode, the
`#ifdef _TIM` AMReX arm with `%to_c`, the `bind(C)` interface, and halo-check/clock relocation.

## Help message

If the arguments are empty, `help`, `--help` or `-h`, print this and stop:

```
Usage: /generate_cpp_bridge <work-directory> <function-name> [--enable_src_validate] [--enable_git_commit] [--disable_git_commit]

Wrap a container-based MOM6 subroutine in a runtime-dispatched shim that
selects between the original Fortran code, a capture mode for offline
validation, and a C++/AMReX bridge.

PRECONDITION: the subroutine's array dummies are already RealArray_t /
IntArray_t containers with a Box_t iteration domain. If not, run
/convert_array_containers first -- this skill does not convert data structures.

Arguments:
  <work-directory>       An existing TURBO-ESM/MOM6 checkout (src/, config_src/),
                         on dev/turbo-debug. Cloning is not performed.
  <function-name>        The subroutine to wrap (case-insensitive).
  --enable_src_validate  Also run Step 1: check the checkout, branch and target.
  --enable_git_commit    Commit and push this run (overrides the global preference).
  --disable_git_commit   Don't commit this run (overrides the global preference).

Without either commit flag, ~/.claude/preferences.json's "git_commit_and_push"
("auto" or "manual"; missing or unreadable means "manual") decides.

Example:
  /generate_cpp_bridge /glade/derecho/scratch/sunjian/MOM6 PPM_limit_pos
```

## Step 0 — validate inputs

Stop on the first failure with a one-line, actionable error; never assume a default or create
anything.
1. Both positional arguments present (else `Error: missing arguments. Run "/generate_cpp_bridge
   --help" for usage.`); the work directory exists.
2. Flags recognized (else `Error: unknown option "<value>"`); `--enable_git_commit` and
   `--disable_git_commit` not both given.
3. `cpp_bridge_lessons` already invoked this session; if not, invoke it now. Never reconstruct a
   template from memory.
4. The target is container-based: every grid-shaped array dummy is `RealArray_t`/`IntArray_t` and
   the iteration domain is a `Box_t`. Otherwise stop with
   `Error: "<function-name>" still takes raw arrays. Run "/convert_array_containers <work-directory> <function-name>" first.`

Items 3–4 always run. No wrapping happens until Step 0 passes.

**Decisions** (take the user's or the campaign plan's value if given; otherwise these defaults):
bridge prefix = the one existing `_bridge) bind(C)` declarations use, else `turbotmp_`
(lessons §3.3); env var = `<UPPERCASE_NAME>_MODE`; only the AMReX arm is gated (lessons §8).
The shim's argument list is not a decision: it is the converted subroutine's, verbatim.

## Procedure

### 1. Validate the checkout (only with `--enable_src_validate`)

- `git remote -v` names `TURBO-ESM/MOM6`; `src/` and `config_src/` exist.
- The tree is clean (never stash or discard). HEAD is at `origin/dev/turbo-debug` (fetch first
  only if the user's instructions allow contacting the remote); if not, ask before switching
  branches.
- Exactly one `subroutine <function-name>` declaration under `src/` or `config_src/`; none is an
  error, more than one means ask which.
- `.claude/skills/generate_cpp_bridge/lessons.md` exists (an install check only; it is a pointer,
  never read it for guidance).
- Print the plan: the file, env var, bridge symbol `<prefix>_<name>_bridge`, and the caller files.

### 2. Map each dummy to its `bind(C)` counterpart

Record intent, kind, rank and optional for every dummy and map it by lessons §3.1–§3.2
(containers → `_C` mirrors, scalars by value with C kinds, pointers → `type(c_ptr)`).

A derived-type dummy with a `bind(C)` mirror in the module (a bundle, shadow or container group)
crosses as that mirror: the module's `<type>_to_c(x)` in the AMReX arm, the mirror type in the
interface (`import ::`). It is captured with the module's `record_<type>` helper: once with no
suffix, then `_before`/`_after` by intent, as for a container (lessons §13). `grid_core_type`
crosses as `Gcore%C` and is captured as `_<field>` for each field the `_fortran` body views,
following `Gcore` into every callee it is passed to, as for a flattened dummy (vertvisc
precedent; a tree root that only forwards `Gcore` still needs its callees' fields). The host-context
`c_ptr` (`convert_calltree` Scope) crosses by value like `OBC` and is captured as
`_host_associated`. A group missing its converter or capture helper is a Phase 2 gap
(`convert_calltree` Step 6b): stop, never write one here.

A derived-type dummy with no mirror is **flattened**: one channel per field. Build the field list
from this subroutine's own body (lessons §10 #5):
1. Grep for every `<dummy>%<field>`, including call actuals; where the whole struct is passed to a
   callee, follow it and repeat. Deduplicate.
2. Classify each field as read and/or written, and record its guards; a field behind its own
   guard (`%associated()` inside the struct's presence check) gets its own null-able channel.
3. Every field gets an interface argument (Step 4), a `%to_c()` in the AMReX arm and capture
   records (Step 5). Print the list with its classification.

### 3. Rename the original

`<name>` → `<name>_fortran` on the `subroutine`/`end subroutine` lines; nothing else in it changes
(lessons §12).

### 4. Add the `bind(C)` interface

At the top of the host module, per lessons §3.3, every dummy doc-commented.

### 5. Write the shim `<name>`

Same name and dummy list as the original, copied verbatim (names, types, intents, `!<` comments),
so call sites need no change. One `select case (mode)` over
`getenv_mode("<ENV_VAR>", default=TIMH_runFORTRAN)` with a capture arm, the AMReX arm in
`#ifdef _TIM` (all `%to_c()` calls inside it), and `case default` calling `<name>_fortran`
(lessons §2, §3.1, §8). Capture names and the `kernel` string follow lessons §13. Casts and null
pointers follow lessons §3.2; check new local names case-insensitively against module names
(lessons §10 #6).

### 6. Update the module imports

Add lessons §14's `use`s, deduplicated.

### 7. Confirm callers need no change

Grep every `call <name>(` and confirm argument count and order still match. Edit no call site; a
mismatch means the shim's signature drifted, so fix the shim.

### 8. Relocate halo checks and clocks

Halo-sufficiency checks move out of the kernel to the caller or shim entry (lessons §15); clocks
stay at the caller (lessons §16).

### 9. Verify

- **Field coverage, by script:** for each flattened dummy, the field set from re-grepping
  `<name>_fortran` and its followed callees equals the interface arguments, the AMReX arm's
  `%to_c()` set, and the capture records (`_after` for every written field). For a mirrored group,
  every field the body uses is in its converter and capture helper, and every written group has an
  `_after` record; for `Gcore`, every field viewed by the body or a callee it passes `Gcore` to is
  recorded. Count/order checks and CAPTURE
  replay can't catch a field missing from all of them.
- **Unchanged code:** `<name>_fortran` differs from the original only by its name lines; no call
  site changed; no existing comment edited.
- **Declarations:** no undeclared or case-insensitively clashing names; `#if` depth balanced.
- **Build and run:** don't compile, build or run tests unless the user explicitly allows it in
  this conversation. Hand the user the lessons §17 matrix to run and say it is unrun. If only the
  shim and capture are delivered, the C++ side of the bridge is the next deliverable.

### 10. Commit (gated)

`--disable_git_commit` skips; `--enable_git_commit` runs; otherwise `git_commit_and_push` in
`~/.claude/preferences.json` decides (`"auto"` runs; anything else, or unreadable, skips and lists
the modified files). Commit everything on `claude_<lowercased name>_bridge` (`git checkout -B`),
with a message naming the subroutine, the bridge symbol, the files touched and the verification
reached, ending with the attribution trailer; push with `-u`. If the push is rejected, stop and
report; never force-push.

## Invariants (hard rules)

- The public name and argument list never change; default-mode numerics are bit-identical.
- No container conversion work (see top).
- Per-kernel env vars only; capture is never gated by `_TIM`; `case default` is always present.
- Never re-implement `turbotmp_helperF` helpers; never reuse a `kernel` string.
- Flattened field lists come from this body, never from another bridge, and are checked by script.
- Doc comments are copied verbatim; no existing Fortran comment is reworded.
- Never claim an unrun build passed.
- Exactly one `!!SKILLS: 0.4.0` marker per touched Fortran file, on its own line after the license
  header and before `module` (update in place). It is the family's shared version.
- Anything not covered: check lessons §10, then ask; never improvise.

## Output to the user

1. The env var and its three values; the `kernel` string.
2. The bridge symbol(s) the C++ side must implement, with any flattened field list.
3. Files touched, and confirmation that no call site changed.
4. The Step 9 results, build status, and what remains on the C++ side; committed, or the files to
   commit.
