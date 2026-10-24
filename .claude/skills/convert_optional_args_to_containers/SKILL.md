---
name: convert_optional_args_to_containers
version: "0.4.0"
description: Convert only the named optional array arguments of a MOM6 subroutine in the middle of a call tree (the pivot) to RealArray_t/IntArray_t containers, leaving its mandatory arguments and its still-raw children untouched. A guarded %view gives a disassociated pointer that is forwarded unconditionally to the child's raw optional dummy, so no if(present(...)) branch trees are needed. Use when converting an optional argument bottom-up would force combinatorial branching at the caller, or when a full convert_array_containers pass isn't wanted yet.
user-invocable: true
argument-hint: <work-directory> <function-name> <optional-dummy-name>[,<optional-dummy-name>...] [--enable_git_commit] [--disable_git_commit]
---

# Convert a subroutine's optional arguments to containers, without converting its children

## Why this exists

If a child's optional dummy becomes a container while its caller (the pivot) stays raw, the
pivot must build a fresh container from its own raw optional dummy; that local is always
present, so "absent forwarded to optional stays absent" no longer holds and the pivot has to
enumerate every present/absent combination (a real case produced a 16-way tree per call site).
Instead, convert the **pivot's** optional dummies and thread containers through its callers;
the pivot forwards into its raw children with no branching. Converting the children is a
separate `convert_array_containers` run.

Campaign decisions recorded in the calling plan (`convert_calltree`) are applied, not re-asked.
Assumes `array_container_lessons` (§6, §9 #8b, §9 #16) is loaded this session.

## Help message

If the arguments are empty, `help`, `--help` or `-h`, print this and stop:

```
Usage: /convert_optional_args_to_containers <work-directory> <function-name> <optional-dummy-name>[,<optional-dummy-name>...] [--enable_git_commit] [--disable_git_commit]

Convert the named optional array dummies of one subroutine to containers,
forwarding them to still-raw children through guarded %view pointers.
Mandatory dummies and children are left untouched.
```

## The mechanism (read this before doing anything)

### Direction 1 — forwarding OUT of the pivot, into a still-raw child

```fortran
real, dimension(:,:), contiguous, pointer :: uhbt   ! named like the child's raw dummy
nullify(uhbt)
if (present(uhbt_a)) call uhbt_a%view(uhbt)
call still_raw_child(..., uhbt=uhbt, ...)
```

A disassociated pointer passed to an `optional` dummy that is neither `pointer` nor
`allocatable` reads as not present. Two conditions make it hold:
- **`nullify` as an executable statement**, before the first guard that might view the pointer.
  A local pointer's status on entry is undefined, so without it the child may see garbage as
  present (a confirmed silent-corruption bug; lessons §9 #16). `=> null()` in the declaration
  is wrong too: it implies `SAVE`.
- **The child's dummy is a plain optional array.** A `pointer`/`allocatable` dummy sees the
  pointer as present-but-disassociated. Check each child's declaration.

### Direction 2 — forwarding INTO the pivot, from its external callers

Decide per call site and per argument, from the caller's actual:
- **A fixed source** (a plain array or expression, or the keyword simply omitted): no branching.
  Build the container unconditionally with `source=`, pass by keyword, `copy2F` if written, `free`.
- **The caller's own optional dummy:** the ordinary guarded-`alloc` problem (lessons §6, §9 #8b);
  a fresh container has no "absent" form. Three or more such optionals at one call site means
  combinatorial branching: stop and ask.

Getting this wrong is caught only by a test that exercises the absent path.

## Steps

### 1. Scope the conversion

Record the pivot and which optional dummies convert; its mandatory and derived-type dummies stay
as they are. List every child each named argument reaches, checking the Direction 1 precondition
(stop on a failure), and every external caller (`call <pivot>(`, including continuations).

### 2. Convert the named optional dummies on the pivot's signature

Intents as in lessons §6: written → `type(RealArray_t), optional, intent(inout)` (never
`intent(out)`), `intent(in)` stays. Rename with `_a`; keep every `!<` comment verbatim. Rename
every `present(<name>)` in the body, grouped tests included; a group partner outside this pass's
scope means the scope is wrong: stop and reconsider.

### 3. Rewire internal forwarding into still-raw children (Direction 1)

One pointer local per converted argument, `nullify`d together before anything else runs; a
guarded `%view` and a keyword forward at every internal call site (all branches, both
directions). Children are not edited.

### 4. Update every external call site (Direction 2)

Apply Direction 2 per site and argument. Arguments outside this pass pass through unchanged.

### 5. Verify

Scripted:
- Every child byte-identical to before.
- Every Direction-1 pointer nullified before any guard that views it.
- Argument count and keywords at every internal and external call.
- `!<` comments preserved verbatim; no loop-body edit; no missed call site; line limit (lessons §7).
- **Build:** don't compile, build or run tests unless the user explicitly allows it in this
  conversation. The mechanism rests on a standard rule inspection can't confirm: say the build
  is unrun and name a test of both the present and absent paths as the first thing to run.

## Versioning marker

Exactly one `!!SKILLS: 0.4.0` line per touched Fortran file, after the license header, before
`module`; update in place if present.

## Hard rules

- No child signature or body edited; no Direction 1 into a `pointer`/`allocatable` dummy.
- Every Direction-1 pointer `nullify`d as an executable statement, never `=> null()`.
- A caller's own optional dummy is never built into an unconditional container.
- Mandatory dummies are `convert_array_containers`'s job; `present()` → `%associated()` is
  `convert_present_to_associated`'s, afterwards.
- No comment about the conversion; a genuinely new dummy still gets a physical-meaning `!<`.
- Never claim an unrun build passed.

## Commit gating

As `convert_array_containers` Step 8. Branch
`claude_<lowercased function-name>_optional_containers`.

## Output to the user on success

Converted dummies (and any left out, with why); internal call sites and each child's precondition
check; external call sites and which Direction 2 case applied; children-unchanged confirmation;
build status with the unverified-mechanism callout; committed, or the files to commit.
