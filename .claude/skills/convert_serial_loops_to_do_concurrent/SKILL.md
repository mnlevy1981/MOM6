---
name: convert_serial_loops_to_do_concurrent
version: "0.4.0"
description: Convert the serial `do` loop nests of MOM6 subroutines to `do concurrent` where a scripted dependency check proves every iteration independent -- pointwise maps become one `do concurrent`, per-point reductions over the outer index are interchanged (`do concurrent` over the points, the reduction serial inside) -- leaving anything unproven serial. Loop bodies never change, so results are bit-for-bit. Use before convert_loops_to_box_iterators on code nobody has yet converted, so every computational loop can map onto an AMReX `ParallelFor`.
user-invocable: true
argument-hint: <work-directory> <function-name>[,<function-name>...] [--enable_git_commit] [--disable_git_commit]
---

# Convert serial loop nests to `do concurrent`

`convert_loops_to_box_iterators` rewrites `do concurrent` headers only; a serial `do` nest is out
of its scope, and an AMReX port needs every computational loop as a `ParallelFor`. This skill
converts a serial nest only when a script proves its iterations independent. It changes loop
control statements only: the body, its order and its arithmetic are untouched, so the result is
bit-for-bit on CPU. What CPU testing cannot show is a dependency or a shared temporary that
would race on a GPU or under OpenMP; that is why the proof is a script, not a test.

Campaign decisions recorded in the calling plan (`convert_calltree`) are applied, not re-asked.

## Help message

If the arguments are empty, `help`, `--help` or `-h`, print this and stop:

```
Usage: /convert_serial_loops_to_do_concurrent <work-directory> <function-name>[,<function-name>...] [--enable_git_commit] [--disable_git_commit]

Convert serial do-loop nests whose iterations are provably independent to
do concurrent (maps) or interchange them (per-point reductions). Bodies untouched.
```

## Classification (`check_loop_independence.py` in this directory)

Run `check_loop_independence.py <file> <subroutines>`; it classifies every serial nest (`do`
statements directly nested, on one line or consecutive lines) and is conservative: anything it
cannot prove is SEQUENTIAL.

| Class | Condition (all must hold) | Conversion |
|---|---|---|
| **MAP** | every array written is indexed by exactly the loop indices (or a scalar that is an index shift, `k = kstart + kk - 1`), with `:` allowed for whole dimensions; an array that is written is read only at the written point; every scalar assigned in the body is assigned before it is read; no `exit`, `return`, `goto`, I/O, or `call` (a `pure`/`elemental` call is confirmed by hand) | one `do concurrent` over the nest's indices, in the original nest order; the assigned scalars go in `DO_LOCALITY(local(...))` |
| **REDUCTION** | as MAP, except some targets omit the **outermost** index and are only accumulated (`X(i,j) = X(i,j) + …`) | interchange: `do concurrent` over the inner indices, the outer index as a serial `do` inside. Each point still accumulates in the original order, so the sum is bit-for-bit. |
| **SEQUENTIAL** | anything else (carried scalars, shifted reads of a written array, inner loops, calls, a write to a scalar derived-type component, an `!$omp`/`!$acc` directive on or in the nest) | none; inner nests of a sequential loop are classified on their own (e.g. a smoothing pass loop stays serial, its sweeps convert) |

Raw-array code works the same: derived-type designators (`CS%Kh_bg_xy(I,J)`, `G%IdxCu(I,j)`) are
treated as arrays, and a nest carrying a device directive stays serial until the plan's
device-directive policy says how to treat it (converting it changes what runs on the device and
what data must be there). `cycle` is allowed in `do concurrent`; `exit` is not. Distinct array names are assumed not to
alias; that holds for `%view` pointers of distinct containers and must be confirmed otherwise.
**Show the table before editing.** A nest the plan expects converted that comes out SEQUENTIAL is
a question for the user, never forced.

## Steps

### 0. Validate inputs

As `convert_array_containers` Step 0 items 1–4; invoke `array_container_lessons` if not loaded.
The file must already include `do_concurrent_compat.h` (for `DO_LOCALITY`) if any converted nest
needs a `local` list; otherwise add the `#include` beside the module's other includes.

### 1. Classify and show

Run the checker on every target; table subroutine, line, nest, class, `local` list or reasons.

### 2. Rewrite the loop control only

- **MAP:** the `do` statements become one header,
  `do concurrent (kk=1:kmax, j=js:je, i=is:ie) DO_LOCALITY(local(k))` (bounds `a,b` → `a:b`,
  index names, case and order unchanged), and the matching `enddo`s become one `enddo`. Text that
  shared a line with the loop statements (`… ; if (mask(i,j) > 0.0) then`, `endif ; enddo ;
  enddo`) keeps its statements; only line breaks around it may change. Never put a statement
  after `DO_LOCALITY(...)` on the same line: the macro expands to `;` on compilers without
  locality support.
- **REDUCTION:** `do concurrent` over the inner indices (with the `local` list), then the outer
  `do` (unchanged) inside it, the body, its `enddo`, and the outer `enddo`.
- A loop stride (`do s=1,0,-1`) or a non-unit step is never converted.

### 3. Verify, by script

- **Bodies:** with loop-control statements removed, every subroutine's statement sequence is
  identical to before (canonicalised for `;`/line breaks).
- **Re-check:** on the edited code the checker reports no MAP or REDUCTION nest; its SEQUENTIAL
  set is the original one plus the serial inner loop of each REDUCTION (serial by design).
- **Locality:** every scalar assigned inside a converted loop is in its `local` list.
- **Format:** line limit (lessons §7); `#if` depth balanced.
- **Build:** don't compile, build or run tests unless the user explicitly allows it in this
  conversation. Say it is unrun; suggest bit-for-bit runs reaching each converted nest. CI does
  not cover GPU offload: say so. Concurrency itself can be exercised on CPU (e.g. nvfortran
  `-stdpar=multicore`, gfortran `-ftree-parallelize-loops=N`): run-to-run differences there point
  at a missed dependency or `local`.

## Invariants (hard rules)

- Only loop-control statements change; bodies, their order and arithmetic never do.
- No nest converted without the checker's MAP/REDUCTION verdict; no SEQUENTIAL nest forced.
- Every scalar assigned in a converted loop is `local`.
- No conversion comment added; never claim an unrun build passed.
- Exactly one `!!SKILLS: 0.4.0` marker per touched file, after the license header.

## Commit (gated)

As `convert_array_containers` Step 8, on branch
`claude_<lowercased first function-name>_do_concurrent`.

## Output to the user

The Step 1 table; per subroutine, nests converted (MAP/REDUCTION) and left (with reasons);
Step 3 results; build status; committed, or the files to commit; `convert_loops_to_box_iterators`
as the next step.
