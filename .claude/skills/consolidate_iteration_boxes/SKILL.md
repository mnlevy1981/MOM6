---
name: consolidate_iteration_boxes
version: "0.4.0"
description: Phase 2 cleanup of a MOM6 call tree's iteration boxes, after convert_loops_to_box_iterators and before bridging. The tree root builds a handful of base boxes (one per stagger, plus any halo decision other kernels depend on, plus any space shared by two or more kernels); every other loop box is one commented grow of a box the kernel receives, named for the points it covers and what it computes. Loop ranges never change, so results are bit-for-bit. Use when kernels build many offset-named boxes (bxU_mi_m2j_pi_p2j) or grow boxes of one stagger to cover another.
user-invocable: true
argument-hint: <work-directory> <root-function> --box <name> | --kernels <k1,k2,...> [--enable_git_commit] [--disable_git_commit]
---

# Consolidate a call tree's iteration boxes

After `convert_loops_to_box_iterators`, every loop iterates a box, but each kernel builds its own
boxes and names them by how they were grown (`bxH_m2ij_pi_p2j`). The halo each loop needs is then
decided deep in the tree, far from the code that decides how much valid halo there is, and the
names hide which points a loop covers. This skill puts the iteration spaces back in order without
changing a single loop range. The model is continuity_PPM: the root grows its base box once by a
stencil width it decides; a few descendants make one commented grow; most use what they receive.

Campaign decisions recorded in the calling plan (`convert_calltree`) are applied, not re-asked.
Assumes `array_container_lessons` (§5 `Box_t` API and memory modes) is loaded this session.

## Help message

If the arguments are empty, `help`, `--help` or `-h`, print this and stop:

```
Usage: /consolidate_iteration_boxes <work-directory> <root-function> --box <name> | --kernels <k1,k2,...> [--enable_git_commit] [--disable_git_commit]

Reorganise the iteration boxes of the call tree under <root-function>: base boxes at the root,
one commented grow per local box, purpose names. Loop ranges are unchanged (bit-for-bit).

  --box bxQstr       Introduce one root box and replace its local grows in every kernel.
  --kernels a,b      Clean up the local boxes of these kernels (names, comments, Rule 4).
```

## Rules

1. **Never change a loop range.** Every loop iterates exactly the points it did before, in both
   memory modes and every option setting. Widening, merging or unifying ranges (including making
   two branches that fill the same array agree) is a separate step with its own proof that the
   computational domain is unchanged; a range that turns out to be a bug is fixed in its own PR.
2. **The root builds the base boxes:**
   - one per stagger (`bxH`, `bxQ`, `bxU`, `bxV`), from the grid's index ranges;
   - any box that encodes a halo decision several kernels must agree on (e.g. how far out the
     viscosity coefficients are computed), built as **one grow of a base box** by widths the root
     computes from the options and the memory mode (`is - Isq` is 1 in symmetric memory and 0 in
     non-symmetric), each width commented with why;
   - any iteration space used by two or more kernels.

   Every root box has a comment giving its range, why it has that halo, and the kernels using it.
3. **Everything else is local** to the one kernel that uses it: **one** `grow`/`growLo`/
   `growHi`/`growBy` call on a box the kernel receives, with constant offsets or widths derived
   from received boxes or control-structure flags, and a comment giving the range, the field(s)
   computed, and why it is that wide. Built once per subroutine and freed on every exit path.
4. **Check-in:** a subroutine that would keep **more than one** local box is shown to the user
   (each box's range, fields, reason and proposed name) before it is edited. There is no hard limit;
   the user decides which stay local and which move to the root.
5. **Names:** `bx` + the stagger of the points the loop computes (`H`, `U`, `V`, `Q`) + a short
   purpose word: `bxUdel2`, `bxQstr`. Base boxes keep their names. Take the stagger from the
   loop's range, not from how the written array is allocated: a scratch array may hold h-point
   values in one kernel and q-point values in another.
6. **Base a box on its own stagger** when that reproduces the range with a constant grow in both
   memory modes. When it doesn't (the original loop is written in another stagger's indices, so the
   width would differ by one between memory modes), keep the base that does reproduce it, and say
   so in the comment ("u points; built from bxH because the range starts at is-2").
7. **Exceptions:** boxes built per pass inside a Fortran orchestrator loop that is never bridged
   (a smoothing pass count that sets the halo) stay where they are.

## Steps

### 0. Validate inputs

As `convert_array_containers` Step 0 items 1–4. The tree must already iterate boxes everywhere
(`audit_bridge_ready.py` loops PASS); otherwise run `convert_loops_to_box_iterators` first.

### 1. Inventory, by script

For every loop in scope (`--box`: every loop whose range equals that box's, in every kernel;
`--kernels`: every loop of those kernels): kernel, line, box, the box's range in both memory modes and each option setting that
changes a root box, the arrays it writes, and how many kernels use the same range. For each range,
the grow of its own stagger box that reproduces it, or "differs between memory modes". **Show the
table before editing**, with the proposed root/local split and names; apply the Rule 4 check-ins.

### 2. Name-clash check

Every new name against every symbol in the touched kernels and the module, ignoring case
(`array_container_lessons` §9 #18). On a clash, stop and ask.

### 3. Edit

Root boxes: declared, built and freed in the root beside the existing ones, passed as `Box_t`
arguments to the kernels that use them (leading the argument list in `bxH, bxQ, bxU, bxV` order,
then the others). Local boxes: declared after `! Local variables`, built once after the
extraction lines, freed before `end subroutine` and every early `return`. Replace only the box
name in each loop header; drop boxes and box dummies left unused. A script may apply this, matching
each header by its exact Step 1 text.

If a kernel that is already bridged gains or loses a box argument, its dispatcher, capture records
and `bind(C)` interface change in the same step (`generate_cpp_bridge`); a rename inside a kernel
needs none of that.

### 4. Verify, by script

- **Ranges:** every loop's range, evaluated from the new boxes, equals the old one in both memory
  modes and every option setting that changes a root box.
- **Bodies:** loop bodies and every statement other than box declarations, builds, frees and the
  box name in loop headers are unchanged.
- **Boxes:** every local box is one grow call on a received box and has a comment; every root box
  has a comment; no box is built twice in one subroutine; every built box is freed on every exit.
- **Calls:** argument count, order and kind match at every call site (and in the bridge triple, if
  touched).
- **Format:** line limit and no trailing whitespace (lessons §7); `#if` depth balanced;
  `audit_bridge_ready.py` still prints BRIDGE-READY.
- **Build:** don't compile, build or run tests unless the user explicitly allows it in this
  conversation. Say it is unrun; ask for bit-for-bit runs in **both** memory modes, with the
  option settings that change the root boxes.

## Staging

One checkpoint per step, each bit-for-bit in both memory modes:
1. **Root boxes first, one `--box` step each.** A step builds the box in the root, passes it to
   every kernel whose loops use that range, removes their local grows of it, and re-bridges every
   kernel whose signature changed (`generate_cpp_bridge`, "Changing the signature of a bridged
   kernel"). This is where kernels lose their box arithmetic.
2. **Then `--kernels` steps** for what stays local: names, comments, the Rule 4 check-ins.

Don't stage by the stagger of the arrays written: scratch arrays reused at two staggers make
that grouping ambiguous.

## Hard rules

- No loop range changes; no loop body changes.
- Every local box: one grow call on a received box, commented. More than one per subroutine is a
  check-in with the user, never decided alone.
- Root boxes built only at the root; no box arithmetic outside the root and the local grows,
  except Rule 7.
- No conversion comments ("renamed from …"); comments say what the box is for.
- Exactly one `!!SKILLS: 0.4.0` marker per touched Fortran file, after the license header.

## Commit gating

As `convert_array_containers` Step 8. Branch `claude_<lowercased root-function>_boxes`.

## Output to the user

The Step 1 table with the root/local split and names; the check-ins and their decisions; the
Step 4 results and build status; committed, or the files to commit.
