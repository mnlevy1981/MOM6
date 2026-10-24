---
name: convert_present_to_associated
version: "0.4.0"
description: Remove `optional` dummies for bind(C) bridge-readiness (C has nullable pointers, not optional arguments). An optional container dummy (RealArray_t/IntArray_t/LogicalArray_t) or container group (shadow/bundle type) becomes non-optional and is tested with `%associated()`; an optional scalar becomes mandatory, with every caller that omitted it passing the value absence used to mean. The unit is an argument name across its whole by-value forwarding chain, not one subroutine. An array argument must already be a container everywhere in that chain (run convert_array_containers first). Adds no bind(C) bridge -- that is generate_cpp_bridge's job.
user-invocable: true
argument-hint: <work-directory> <dummy-argument-name>[,<dummy-argument-name>...] [<seed-function-name>]
---

# Convert `present(x)` to `x%associated()`

## Why this exists

`bind(C)` can't express an omitted argument. The argument becomes **always present as a
container** whose data may be unassociated, tested with `%associated()` (never `x%data`, which
leaks the field name out of `array_mod`). This is for bridge-readiness only; both tests are O(1).

Campaign decisions recorded in the calling plan (`convert_calltree`) are applied, not re-asked.

## Help message

If the arguments are empty, `help`, `--help` or `-h`, print this and stop:

```
Usage: /convert_present_to_associated <work-directory> <dummy-argument-name>[,<dummy-argument-name>...] [<seed-function-name>]

Drop `optional` from a dummy across its whole forwarding chain. A container
dummy is tested with x%associated() and callers that omitted it pass an
unallocated placeholder; a scalar dummy becomes mandatory and callers that
omitted it pass its absent-default explicitly.
```

## Step 0 — scope: one subroutine, or the whole forwarding chain

The presence idiom must match at both ends of every **by-value** forward (a guarded `%view`
carries its own signal and ends the chain). Build the set by closure:

1. **Seed** from `<seed-function-name>` (confirm it declares the argument as an optional
   container), or from every subroutine declaring an optional dummy of that name that is the
   same quantity.
2. **Close** it, repeating until a pass adds nothing, from *every* member (not just the seed):
   *down* — a callee whose same-named optional dummy receives it by value joins; *up* — a caller
   that passes its own same-named optional dummy joins. A chain ends where the argument is
   consumed (`%view`/`%associated()`) or passed to a mandatory dummy under another name. A
   caller that passes a fixed value or omits it is a **boundary**: not a member, but its call
   site gets Step 4.
3. Members get Steps 1–3; every boundary call site gets Step 4, in the same pass. A one-member
   set is common, but only after the closure says so.

## Hard precondition

An array argument is already a container in every member (a scalar or derived type needs nothing). A `present()`-linked group (`present(a) .and.
present(b)`, or all-or-nothing handling) converts together, never one member alone.

## Steps

1. **List every `present(<name>)`** in every member before editing.
2. **Drop `optional`** from each declaration; nothing else in it changes.
3. **Replace `present(x)` with `x%associated()`**, checking each occurrence really tests this
   argument. `%associated()` must be `pure` (it may sit inside `do concurrent`); confirm once
   per container type. At every rewritten guard, independently:
   - a pointer local set by a guarded `%view` is `nullify`d before the first guard (lessons §9 #16);
   - inside an `!$omp target`/`teams`/`parallel`/`loop`-family construct, `%associated()` is never
     referenced: precompute a logical before the construct (lessons §9 #17).
4. **Fix every boundary call site.** One that already passes a container needs nothing. One that
   omitted the argument passes a never-`%alloc`'d placeholder local, declared inside the
   subroutine that contains the call (check its line range, not proximity). One placeholder never
   serves two dummies of one call when either is `intent(out)`/`inout`: that aliasing is illegal
   whichever branch runs.
5. **Leave `%to_c()` alone**; it is already null-safe.

## Optional scalar dummies

An optional `real`/`integer`/`logical` scalar can't cross `bind(C)` either, and has no
"unassociated" state, so it becomes **mandatory** and its absent meaning moves to the callers.
Step 0's closure applies unchanged (a caller forwarding its own same-named optional joins the set).

1. **Find the absent-default:** the single value the body uses when `present(x)` is false,
   typically `x_val = <default> ; if (present(x)) x_val = x`. If absence selects behaviour no
   single value expresses, stop and ask.
2. **Drop `optional`.** Replace only the fallback: the default-then-override lines become
   `x_val = x`; the rest of the body is unchanged.
3. **Every boundary call site that omitted it passes the default explicitly** (by keyword if the
   call already uses keywords); call sites that passed a value keep it.
4. An existing doc comment that describes the absent case ("If present and false, …") is left
   verbatim and reported to the user (lessons, comment STANDING RULE).

## Optional derived-type dummies (shadows, bundles)

A derived type grouping containers and scalars (a `create_shadow_container_type` shadow, a
bundle) is handled like a container: **drop `optional`**; the caller always passes an instance,
and a caller that omitted it passes a never-built instance (its containers unassociated, its
scalars at their default component values). Each `present(x)` guard becomes the presence test of
what it protects:
- around a container field's use: that field's `%associated()`, or nothing where an inner guard
  already tests the field (keep the `nullify`, lessons §9 #16);
- around a scalar field's read: no test, provided the type's default component value equals the
  old absent-default; otherwise stop and ask. A default-then-override pair becomes the plain read.

## Verify

Scripted:
- No `present(<name>)` or `optional … <name>` left repo-wide for any argument in the set.
- Scalars: every call that omitted the argument now passes the recorded absent-default.
- Every placeholder declared within its call's enclosing subroutine; no actual repeated in one call.
- Every guarded-`%view` pointer nullified first; no `%associated()` inside an OpenMP construct.
- Doc comments verbatim; line limit (lessons §7).
- **Build:** don't compile, build or run tests unless the user explicitly allows it in this
  conversation; say it is unrun.

## Explicitly not performance motivated

`present()` and `%associated()` are both O(1) runtime checks; the only gain is bridge-readiness.

## Versioning marker

Exactly one `!!SKILLS: 0.4.0` line per touched Fortran file, after the license header, before
`module`; update in place if present.

## Hard rules

- Never convert a raw array, one member of a linked group, or a subroutine without the closure.
- A scalar's absent-default is one value read from the body, never guessed; otherwise ask.
- A derived type's scalar field read unguarded only if its default component value is the absent-default.
- Callee and callers change in the same pass.
- Placeholders scoped by line range, never shared across two definable dummies of one call.
- `%associated()` is `pure`; nullify and OpenMP hazards checked at every rewritten guard.
- No `bind(C)` interface, dispatcher or capture mode (`generate_cpp_bridge`, afterwards);
  never claim an unrun build passed.
