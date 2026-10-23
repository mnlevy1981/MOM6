---
name: convert_array_containers
version: "0.3.2"
description: Convert a MOM6 Fortran subroutine from raw grid-shaped arrays to RealArray_t/IntArray_t/LogicalArray_t containers in place -- array dummies become containers, the body gets %view pointers so the math is untouched, and every call site gets alloc/copy2F/free. Also fixes local-array sizing still keyed off G/GV macros once a matching dummy is already a container (Step 2a), and promotes any remaining G%/GV%/US%/CS% field reference, array or scalar, anywhere in an already-converted subroutine's body to a plain dummy (Step 2b). Use when pushing the raw-array boundary up the MOM6 call tree. Does NOT rename anything, add a dispatcher shim, capture mode, or a bind(C) bridge -- those belong to generate_cpp_bridge.
user-invocable: true
argument-hint: <work-directory> <function-name> [--enable_src_validate] [--enable_git_commit] [--disable_git_commit]
---

# Convert a MOM6 subroutine to array containers

This skill is the **execution checklist**. Templates, the container API
reference, intent rules, and pitfalls live in `array_container_lessons`
(§1–§10) — invoke it once near the session start, before this skill, and
refer back to its sections from each step below (Step 0 checks it's
loaded). Do not reproduce templates here.

**Data structures only.** No rename, no `*_fortran` variant, no
dispatcher/`getenv_mode`/capture mode/`#ifdef _TIM`/`%to_c`/`bind(C)`
interface — those are `generate_cpp_bridge`'s job, layered on afterwards
(lessons §10). A routine may be converted with no commitment to bridging.

## Help message

If `$ARGUMENTS` is empty, or equals `help`, or equals `--help`, or
equals `-h`, do NOT run any steps. Print the following help message
verbatim and stop:

```
Usage: /convert_array_containers <work-directory> <function-name> [--enable_src_validate] [--enable_git_commit] [--disable_git_commit]

Convert one MOM6 subroutine from raw grid-shaped Fortran arrays to
RealArray_t / IntArray_t / LogicalArray_t containers, in place. The subroutine keeps its
name; its array dummies become containers; its body obtains pointers via
%view so the math is unchanged; and every call site is updated to
alloc / copy2F / free around the call.

Does NOT rename to *_fortran, add a dispatcher shim, capture mode, or a
bind(C) bridge -- use /generate_cpp_bridge for that, afterwards.

Arguments:
  <work-directory>   Absolute path to an existing TURBO-ESM/MOM6 checkout
                     (must contain src/ and config_src/). Cloning is not
                     performed.
  <function-name>    Name of the Fortran subroutine to convert
                     (case-insensitive match against its declaration).
  --enable_src_validate  (optional) Run Step 1: verify the work directory is a
                         TURBO-ESM/MOM6 checkout with a clean tree, and that
                         the subroutine exists. Off by default.
  --enable_git_commit    (optional) Force Step 9 to run this invocation
                         only, overriding the global git_commit_and_push
                         preference in ~/.claude/preferences.json.
                         Mutually exclusive with --disable_git_commit.
  --disable_git_commit   (optional) Force Step 9 to be skipped this
                         invocation only, overriding the global
                         preference. Mutually exclusive with
                         --enable_git_commit.

When neither --enable_git_commit nor --disable_git_commit is passed, Step 9
follows ~/.claude/preferences.json's "git_commit_and_push" key ("auto" or
"manual"; treated as "manual" when the file is missing, the key is absent,
or the JSON fails to parse).

Example:
  /convert_array_containers /glade/derecho/scratch/sunjian/MOM6 zonal_mass_flux
```

## Step 0 — validate inputs

`$0` = work-directory, `$1` = function-name. Run these checks **before any
other step** and stop on the first failure with a one-line, actionable
error. Do not retry, do not assume defaults, do not create anything.

1. **Argument count.** `$0` or `$1` empty → stop:
   `Error: missing arguments. Run "/convert_array_containers --help" for usage.`
2. **Work directory exists.** If `$0` isn't an existing directory → stop:
   `Error: work directory "<value>" does not exist.` (No cloning here —
   it must already be a checkout.)
3. **Parse optional flags** — `--enable_src_validate`,
   `--enable_git_commit`, `--disable_git_commit`. The latter two are
   mutually exclusive; an unrecognised `--` flag stops with a usage error.
4. **Reference material primed.** Assumes `array_container_lessons` was
   already invoked this session (§2–§9 are referenced throughout). If you
   have no memory of its sections, invoke it now — never guess template
   content or API signatures from a similar prior task.
5. **A `<function-name>_fortran` sibling doesn't mean stop.** If
   `grep -irn "^[[:space:]]*subroutine[[:space:]]\+<function-name>_fortran\b" $0/src $0/config_src`
   matches, `<function-name>` is a bridge shim over an already-converted
   worker — convert it exactly like any raw-array subroutine (Steps 2–9),
   with one constraint (lessons §10): never edit inside
   `select case (mode)…end select` except renaming a `CS%x`/`GV%x`
   reference when Step 2 drops that type, and never touch a `bind(C)`
   interface. Either need signals drift into `generate_cpp_bridge`'s
   territory — stop.
6. **Target isn't already converted — or has Step 2a/2b work left.**
   Locate `<function-name>`'s declaration
   (`grep -irn "^[[:space:]]*subroutine[[:space:]]\+<function-name>\b" $0/src $0/config_src`)
   and read its dummies. If every array dummy is already a container,
   this subroutine's own dummy-conversion (Steps 2, 4, 5, 6, 7) is done —
   but check both, since either or both can still apply:
   - **Step 2a candidate:** a local array (any element type), or a raw
     dummy the user chose to leave unconverted in Step 2 (the only way a
     dummy still reaches this step), still sized via
     `SZI_`/`SZJ_`/`SZK_`/`SZIB_`/`SZJB_` and matching a mandatory
     dummy's shape (Step 2a).
   - **Step 2b candidate:** grep the body for every `<type>%<field>`
     occurrence of a derived-type dummy this subroutine still takes but
     isn't converting (`G`/`GV`/`US`/`CS`/other) — including inside a
     call argument, not just the subroutine's own expressions. Any hit
     not already reached through an existing container dummy is a
     candidate, whether used directly, only forwarded, or both (Step 2b).

   Neither applies → stop, report already container-native. Either or
   both apply → proceed with exactly those, one invocation not two; Step
   3 still runs for call-site tracing, Step 9 still verifies whatever
   changed, every other dummy-conversion step stays skipped.

Items 4–6 always run, regardless of flags. Step 1 runs only with
`--enable_src_validate`. Step 9's behavior follows
`--enable_git_commit`/`--disable_git_commit` if passed, else the global
preference (Step 9). Nothing else runs until Step 0 passes.

## Settle these decisions (ask if not obvious from the tree)

1. **Iteration box** — use an existing `type(Box_t)` dummy (commonly
   `bxC`) if there is one. Otherwise the loop bounds derive from `G`/`GV`
   index fields; ask the user whether to add a `Box_t` dummy or keep the
   existing bounds. Never invent one silently.
2. **Argument-list reduction** — whether `G`/`GV`/`US`/`CS` can drop from
   the signature once their members become containers/plain dummies
   (lessons §8). Default: drop only what the body no longer references —
   but that can't be finalized until every callee this subroutine still
   forwards the type to wholesale has settled the same question itself
   (direction note below); a clean-looking body doesn't mean `G` can go
   if it's still passed to a still-forwarding descendant. If a *private*
   CS's flattened fields keep recurring together across leaves, don't
   leave them flattened — run `create_config_bundle_type` afterward to
   cluster them into small bind(C)-ready bundles.
3. **Scope of call-site updates** — default to updating every call site
   from Step 3 in the same change, so the build stays green. Each site
   is handled per its caller's own state (Step 8, Case A/B).

One subroutine per invocation, no recursion — a multi-level conversion
is a sequence of invocations. Two direction questions apply, and they
needn't agree:

- **Dummy conversion (Steps 2/4/5/7/8): top-down preferred** (lessons
  §1) — the marshalling block lands once, at the permanently-raw caller,
  instead of being written and torn down at every level.
- **Dropping `G`/`GV`/`US` (item 2, Step 2b): bottom-up *required***, not
  just preferred. A subroutine can't safely drop `G` while any callee
  still needs it forwarded whole, and that's only knowable once every
  callee has settled it — deciding top-down means guessing or
  re-checking later. Walk leaves-first so each drop decision is made
  exactly once, correctly.

In practice: containerize dummies in whichever order suits marshalling
churn, but defer every `G`/`GV` drop decision and Step 2b promotion
until walking back up from the leaves. Report Step 7's follow-up
candidates so the user can pick the next invocation deliberately, and
say which direction the current one implies for both questions.

If the user already specified any of these, take their values as-is.

## Procedure

Each step is one action with a pointer to the `array_container_lessons`
section that holds the template or rationale.

### 1. Validate the checkout *(only when `--enable_src_validate` is passed; otherwise skip to Step 2)*

   Verify `$0` is a TURBO-ESM/MOM6 checkout (`git -C $0 remote -v` contains
   `TURBO-ESM/MOM6`, else stop: `Error: "<value>" is not a TURBO-ESM/MOM6 checkout.`).

   - Dirty tree (`git -C $0 status --porcelain` non-empty) → stop and
     surface it; never stash or discard.
   - Report the current branch in the plan confirmation below — this
     skill needs no specific base branch, but the user should see which
     one is about to change. Do not switch branches.

   Then, stopping on the first failure with a one-line error:
   - **Looks like MOM6.** `$0/src` or `$0/config_src` missing → stop.
   - **Subroutine present.**
     `grep -irn "^[[:space:]]*subroutine[[:space:]]\+<function-name>\b" $0/src $0/config_src` —
     0 matches → stop; >1 → list candidates and ask which to convert.
     (Step 0 items 5–6 already ruled out an already-bridged/-converted
     target; no need to re-check.)

### 2. Classify every dummy argument
   Read the full subroutine. For each dummy record: name, declared type,
   `intent`, rank, and whether it's a grid-shaped array, a scalar, a
   `logical`, a `pointer`, or a derived type.

   Decide per dummy (lessons §6, §7):
   - Grid-shaped `real`/`integer`/`logical` array → `type(RealArray_t)`/
     `type(IntArray_t)`/`type(LogicalArray_t)`, renamed `<name>_a` — all
     three element types get the same treatment.
   - Scalar / `pointer` (`OBC`) → unchanged. A bare `logical` *scalar*
     also passes through unchanged; only *arrays* containerize.
   - Derived type (`BT_cont_type`, …) → unchanged (lessons §9 #11).

   **Stop and ask about every `optional` array dummy, of any element
   type — never decide silently.** This skill has produced both
   outcomes in practice (`zonal_mass_flux`'s `uhbt`/`visc_rem_u`/
   `u_cor`/`du_cor` were left raw; other conversions containerized
   theirs) and neither is a safe default. Present each by name; ask:
   convert now, or leave raw.

   - **Convert now.** Becomes a container, keeps `optional`; its `%view`
     is guarded by `present(<name>_a)`, and every existing
     `present(<name>)` renames to match (lessons §6, §9 #8a) — note it
     now so Steps 4, 5, 8 all handle it; a missed guard compiles and
     fails only at run time.
   - **Leave as raw.** No change to its declaration or `present()`
     tests; Step 8 simply forwards it as-is.

   Record every answer before Step 3 — no inferred answers.

   Then scan the body for derived-type array references (`G%mask2dT`,
   …) and scalar references (`GV%Angstrom_H`, …). Arrays become new
   container dummies; scalars become plain dummies computed at the call
   site (lessons §8). **Exception: `G%<array>` references** come from
   `grid_core`/`grid_OBC` (see Step 2b), never a new `<f>_a` dummy.

   Print the proposed new argument list — including any Step 2a
   rewrites — and pause for user confirmation before writing anything.

### 2a. Eliminate grid-macro sizing on locals (or a raw-left dummy) that an existing dummy's bounds already describe

   For every local array (`real`/`integer`/`logical`, any rank), and any
   dummy left raw by an explicit Step 2 choice (an `optional` array the
   user chose not to convert — the only way a dummy still reaches this
   step), still sized via `SZI_`/`SZJ_`/`SZK_`/`SZIB_`/`SZJB_`: check
   whether any **mandatory** array dummy of this subroutine (converted
   here or previously) has the identical macro at some matching dimension.

   - **Match dimension-by-dimension, not rank-by-rank.** A rank-2 local
     `SZIB_(G),SZK_(GV)` matches a rank-3 dummy `u_a` whose first
     dimension is `SZIB_(G)` and third is `SZK_(GV)`, ranks differing
     and the middle dimension skipped — record which dummy/dimension
     covers each of the local's own dimensions; different dimensions may
     key off different dummies.
   - **Only mandatory dummies qualify.** An `optional` dummy's container
     can be absent, so its `%lb`/`%ub` can be disassociated — sizing
     another array off it would be undefined behavior exactly when it
     matters most. If the only macro-matching dummy is `optional`, leave
     that dimension untouched and note it in the Step 10 report.
   - **For a dummy converted in a prior pass**, its raw declaration no
     longer exists — trace one call site (Step 3; they must all agree)
     to the raw array supplied as `source=` in that dummy's `%alloc`,
     and read *that* declaration to recover the macro pattern.
   - Rewrite from macro form to the dummy's own bound fields:
     ```fortran
     ! BEFORE
     real, dimension(SZIB_(G),SZK_(GV)) :: uh_aux

     ! AFTER -- u_a is a mandatory dummy of this subroutine, already a container
     real, dimension(u_a%lb(1):u_a%ub(1), u_a%lb(3):u_a%ub(3)) :: uh_aux
     ```
     Legal because `u_a` is a genuine dummy, associated by the caller
     before entry — unlike a `%view`-derived pointer, only associated
     partway through the body and never legal to reference this way
     (Hard rules).
   - **A mandatory `logical` dummy is not this step's job** — it now
     containerizes exactly like `real`/`integer` (Step 2); don't apply
     this bound-rewiring trick as a substitute. The only dummy still
     eligible here is one left raw by an explicit Step 2 choice.
   - A local/dummy with no macro-matching mandatory dummy is unaffected.

   Do this **before** Step 4 rewrites any dummy declaration — the
   "before" macro pattern is only readable pre-rewrite. If only Step 2a
   applies (no Step 2b), it's the entire run — Steps 4–8 have nothing
   left to do since the signature doesn't change. If Step 2b also
   applies, do this step first (it needs the final mandatory-dummy set).

### 2b. Promote remaining derived-type field references, even on an already-converted subroutine

   Step 2a only rewrites a *bound expression*; it never touches a live
   field reference in the body. This step finishes what Step 2 already
   does on first conversion — extracting a derived type's fields into
   plain dummies — for a subroutine already past Steps 2/4/5.

   **The check is one grep, not a judgment call.** For each derived-type
   dummy this subroutine still takes but isn't converting (`G`/`GV`/
   `US`/`CS`/other), grep the body for every `<type>%<field>` occurrence
   — including inside a call argument, not just the subroutine's own
   expressions. Deduplicate by field name. For every distinct field,
   regardless of how it's used:

   - **Array field** (`G%dy_Cu`, …): if a local container already exists
     in this subroutine sourced from that exact field (typically built
     to forward it into an already-converted callee), promote *that*
     local to a dummy and delete its now-redundant `%alloc`/`%free`.
     Otherwise create a new dummy `<field>_a` with a fresh `!<` comment
     (Step 4). Either way, declare/reuse a `%view` pointer under the
     field's bare name (Step 5), replace every literal occurrence, and
     thread the dummy through every call site via Step 8's ordinary Case
     A/B logic — a raw caller builds/frees a container (flag a
     loop-invariant source as a hoisting candidate); an already-converted
     caller passes its matching dummy straight through or needs its own
     Step 2b pass first (report as a follow-up, don't convert it here).
   - **Scalar field** (`GV%H_subroundoff`, …): add a plain scalar dummy
     with a fresh `!<` comment, replace every occurrence with its bare
     name, and at every call site add
     `<field>=<caller's-own-derived-type>%<field>` — no container.

   **Grid data comes from `grid_core`/`grid_OBC`, never per-field promotion.** Every
   `ocean_grid_type` array a planned tree reads has a persistent copy, built once at init, in one
   of two singletons in `src/core/MOM_grid_containers.F90` (lessons §5, §8):
   - `grid_core`: every non-OBC grid array;
   - `grid_OBC`: `OBCmaskCu IdxCu_OBCmask OBCmaskCv IdyCv_OBCmask`.
   
   The rules:
   - A subroutine that reads `G%<f>` takes `type(grid_core_type), intent(in) :: Gcore` (and/or
     `type(grid_OBC_type), intent(in) :: Gobc`) at `G`'s position.
   - It declares `real, dimension(:,:), contiguous, pointer :: <f>`, calls
     `call Gcore%<f>%view(<f>)` next to its other views, and replaces `G%<f>` → `<f>` in the body.
   - The cap fetches `Gcore => grid_core(G)` / `Gobc => grid_OBC(G)` and passes them down.
   - **No per-call `alloc(source=G%<f>)` and no `<f>_a` dummy.**
   - If a field is missing from both structures, stop and ask the user whether to add it. Never
     promote it per field, and never add a third group without asking.
   - Grid *index scalars* (`G%isc`, `G%IscB`, …) are not in either structure. At a tree root they
     come from boxes the cap builds (e.g. `bxH0`/`bxQ0`), not from scalar dummies.
   - New `Gcore`/`Gobc`/pointer declarations get new comments. The removed `G` declaration goes
     whole; no other comment changes (see the standing rule in `array_container_lessons`).

   Do this for **every** distinct field found — one used directly here is
   exactly as eligible as one only forwarded, and finding one says
   nothing about the others. The derived type drops from the signature
   (Settle these decisions #2) only once the grep comes back empty.

   Replacing a literal `<type>%<field>` with a bare name is the one
   deliberate exception to "zero body edits" — every changed line must
   differ from before by *only* that substitution; verify line-by-line
   in Step 9.

   **Run this bottom-up, leaves first** — a subroutine's dependence on
   these types can't be settled until its callees have settled theirs.

   Present the full field list for confirmation before writing, alongside
   whatever Step 2a found. If only this step applies, skip straight here
   from Step 0/1 — Step 3 still runs for its own call-site threading.

### 3. Find every call site
   `grep -rn "call[[:space:]]\+<function-name>(" $0/src $0/config_src`,
   plus continuation-line forms
   (`grep -rn -A2 "call[[:space:]]\+<function-name>[[:space:]]*&"`).
   List every file/line — this is the Step 8 work queue; a missed site
   breaks the build (lessons §9 #9). A site inside a routine you weren't
   asked to touch is expected: the marshalling block lands there (lessons §1).

### 4. Rewrite the dummy declarations
   Apply Step 2's classification to the declaration block. Written
   containers take `intent(inout)`, never `intent(out)` (lessons §9 #6).

   **Every `!<` doc comment must survive verbatim.** The single most
   commonly broken part of a conversion: the original declaration is
   usually two lines with the comment on the continuation line, while
   the converted one is a single line, so collapsing the two silently
   drops the comment:

   ```fortran
   ! BEFORE -- comment lives on the continuation line
     real,  dimension(SZI_(G),SZJ_(G),SZK_(GV)), &
                              intent(in)    :: h_in !< Tracer cell layer thickness [H ~> m or kg m-2].

   ! AFTER -- WRONG, comment dropped during the collapse
     type(RealArray_t),    intent(in)    :: h_in_a

   ! AFTER -- CORRECT, comment carried onto the new single line
     type(RealArray_t),    intent(in)    :: h_in_a     !< Tracer cell layer thickness [H ~> m or kg m-2]
   ```

   Procedure to guarantee this:
   1. Before editing, record each array dummy's doc comment verbatim,
      including unit annotations (`[H ~> m or kg m-2]`, `[nondim]`).
      Scan the whole declaration — the comment may be on the `real,
      dimension(...)` line or any continuation.
   2. Write the new declaration, re-attach the recorded comment exactly
      as worded — never reword, abbreviate, "improve," or shorten it to
      fit a line; wrap instead (item 5).
   3. A **new** dummy with no prior comment (a grid-derived container, a
      scalar lifted out of `CS`) gets a fresh `!<` comment in house
      style — every dummy here carries one, and Doxygen warns otherwise.
   4. Align the `!<` column with neighbours for readability, never at
      the cost of the text.
   5. **Check the resulting line length — MOM6's limit is 100
      characters.** Collapsing a two-line declaration onto one often
      exceeds it, especially with a unit annotation; wrap onto a `!!`
      continuation at a word boundary rather than shortening anything,
      matching the style neighbouring dummies already use. Recheck both
      resulting lines.
   6. **The `subroutine <name>(...)` header line has the same limit —
      check it separately.** New grid-derived dummies (Step 2) lengthen
      it, and it's easy to forget while checking each declaration. If
      the first physical line exceeds 100 characters, move the `&`
      continuation earlier, breaking the argument list at a comma —
      never remove or rename an argument to fit. Recheck every physical
      line, not just the first.

   Comments not on a dummy declaration (the `!>` header, blocks between
   declarations, anything in the body) stay exactly as they are.

   Do **not** rename the subroutine (lessons §9 #1).

### 5. Add the `%view` locals and calls
   Declare one pointer per container, using the **original** array name
   (matching rank): `real, dimension(:,:,:), contiguous, pointer :: h_in, h_W, h_E`.
   Immediately before first use, `call <name>_a%view(<name>)` for each
   (lessons §2, §4.3).

   For an `optional` dummy flagged in Step 2, guard the view and rename
   its `present` tests (lessons §6):
   `if (present(hin_a)) call hin_a%view(hin)`, and every existing
   `present(hin)` becomes `present(hin_a)`.

   Two mandatory, independent checks on every guard introduced here
   (confirmed real CI bugs — lessons §9 #16-17):
   - Every pointer local set only by the guarded `%view` needs a
     `nullify` before the first such guard, same edit.
   - If the guard sits inside an `!$omp target`/`teams`/`parallel`/
     `loop`-family construct, never leave `%associated()` itself as the
     in-loop condition — precompute a scalar outside the construct.

   The loop body needs **zero edits** apart from those `present`
   renames — the pointers carry the container's original `lb:ub`. Any
   edit to an index expression or formula means something's wrong.

### 6. Wire up the iteration box
   Use an existing `Box_t` as-is. If the body needs a grown/shrunk
   domain, derive it with `bx = bxC%grow(dim=N, n=M)` — `dim=1` zonal/i,
   `dim=2` meridional/j (lessons §5) — and `call bx%free()` when done;
   `grow` returns a new, owned box (lessons §9 #3). Leave existing
   loop-bound expressions alone unless the user chose to add a `Box_t`
   ("Settle these decisions" #1).

### 7. Handle `elemental` callees and unconverted descendants
   An `elemental` scalar kernel called inside a loop (`flux_elem`,
   `flux_elem_OBC`) keeps its signature; pass `%view` pointer elements
   to it exactly as before (lessons §9 #10).

   A callee still taking raw arrays is below the boundary and out of
   scope here: pass it the `%view` pointer, not the container, and note
   it in the Step 10 report as a follow-up candidate.

   If the argument forwarded to that still-raw callee is itself one of
   *this* subroutine's `optional` container dummies, a plain `%view`
   isn't enough — an unguarded `%view` on an absent container is
   invalid, and the callee's `present()` check needs a real signal, not
   just a pointer. Use the guarded-view/mandatory-`nullify`/
   unconditional-forward technique instead (lessons §9 #16); see
   `convert_optional_args_to_containers` for the full mechanism.

### 8. Update every call site from Step 3
   **Classify each call site first.** Look at the *calling* subroutine's
   own dummy declarations and decide which case applies — don't assume;
   a campaign hits both, often in the same file.

   **Case A — the caller still takes raw arrays.** Insert the
   marshalling block (lessons §3): a container per array argument,
   `alloc` every array **with** `source=` — including `intent(out)`
   ones, unless you've confirmed the callee's write covers the full
   `LBOUND`/`UBOUND` rather than just a `Box_t`-scoped sub-region
   (lessons §6) — call the converted routine, `copy2F` outputs/inouts
   back, `free` every container. Derive new scalar dummies here
   (`h_min = 2.0 * GV%Angstrom_H`).

   Common when converting bottom-up; also the permanent end state when
   the caller will never be converted (a module entry point like
   `continuity_PPM`).

   - *Array from the caller's own local scratch variable* (not a
     dummy): this block is itself a `convert_locals_to_containers`
     candidate later — not this skill's job now; note as a follow-up
     once this callee's conversion has settled.
   - *Caller's source array is `optional`*: a container local is always
     present as an actual argument, so an unallocated one makes
     `present()` true in the callee with garbage behind it. Allocate
     conditionally, branch the call (lessons §6, §9 #8b); more than two
     optional arguments at one site multiplies branches — stop and ask
     rather than nesting conditionals.
   - *Loop-invariant source* (grid metadata like `G%IareaT`,
     `G%mask2dT`): do **not** restructure the caller to hoist it as part
     of this conversion. Emit the straightforward per-site `alloc`/
     `free` and **report it in Step 10** as a hoisting candidate — safer
     to do once, with every call site visible, than incrementally.

   **Case B — the caller is already container-based.** No marshalling
   block; pass containers straight through:

   - Matching container already a caller dummy → pass it directly, no
     `alloc`/`copy2F`/`free`.
   - Caller was passing a `%view` pointer here because it used to be raw
     → replace with the container it came from; delete the pointer
     local if nothing else uses it.
   - Caller holds a marshalling block for *this* call from an earlier
     bottom-up pass → now dead scaffolding (lessons §3): delete it, pass
     the caller's own containers.
   - Only a genuinely new array with no existing container (a local
     work array) still needs its own `alloc`/`free`.

   Common when converting top-down — the marshalling is written once at
   the permanent boundary instead of added and removed at every level.

   **Both cases:** don't move existing `cpu_clock_begin`/`cpu_clock_end`
   boundaries; ensure early returns still reach any `free` block
   (lessons §9 #15). Re-check argument count/order at every site against
   the new signature.

   **Raw-vs-container check, argument by argument.** Re-read Step 2's
   list of dummies that stayed raw for *this* subroutine and confirm the
   actual argument at each such position is still raw, not a container
   — easiest to get wrong on a subroutine that mixes both, since a call
   site is usually built by mechanically appending `_a` across most of
   the list, and a raw dummy's caller-side actual often has a same-stem
   container sibling sitting right next to it in scope (a raw `uh_3d`
   fed by the caller's own `uh`, while the caller also holds `uh_a` for
   a different argument). Passing `uh_a` there is a real type mismatch
   only a compiler catches — which this skill can't assume is available
   (Step 9). Caused a real build failure in an already-reviewed
   conversion (`zonal_flux_adjust`'s `uh_3d`), caught only by CI. Check
   this explicitly rather than trusting the mechanical edit.

### 9. Verify
   Run the checks in lessons §10:
   - **Doc comments intact.** Count `!<` before/after per edited
     subroutine; added count ≥ removed (greater if Step 4 introduced new
     dummies). `git -C $0 diff -U0 -- <file> | grep '^-' | grep -c '!<'`
     vs. the `'^+'` version. Confirm each removed comment's text
     reappears somewhere added.
   - **Line length ≤ 100** (Step 4 items 5–6) on every added line,
     including the header:
     `git -C $0 diff -U0 -- <file> | grep '^+' | grep -v '^+++' | awk
     '{ print length($0)-1 }' | sort -rn | head -1`. Wrap, never shorten.
   - Diff review: no edits inside a loop body; name unchanged; only
     declarations, `%view` calls, and call-site blocks changed.
   - Every Step 3 call site updated; argument count/order matches.
   - **Raw-vs-container argument audit** — walk each call site's actuals
     position by position against the callee's current signature;
     matching count/order alone doesn't catch a same-position type swap
     (Step 8).
   - Nothing outside the FMS2 ∩ TIM shared API surface (lessons §4) — no
     `%to_c` (lessons §9 #12).
   - **Guarded-`%view` nullify audit (lessons §9 #16).** Grep for every
     `pointer` local on the left of a guarded `%view` (inside an
     `if (present(...))`/`if (...%associated())` block, not
     unconditional) and confirm a `nullify` sits before the first such
     guard. Add any missing one now, don't defer it.
   - **OpenMP-construct audit (lessons §9 #17).** No `%associated()`
     check, new or pre-existing, is or sits inside an `!$omp
     target`/`teams`/`parallel`/`loop`-family construct's condition —
     every such guard precomputes a scalar outside it. Fix now if not.
   - **Step 2a rewrites, if any:** confirm every rewritten bound keys off
     a genuinely mandatory dummy (recheck even one already a container
     before this run), and references the container itself, never a
     `%view` pointer.
   - **Step 2b promotions, if any:** re-run the field grep — should come
     back empty; a hit means a field was missed. Confirm every changed
     body line differs from "before" by only the substitution; a merged
     array field's old `%alloc`/`%free` is gone; a promoted scalar
     introduced no container/`%view`/`alloc`/`free`.

   **Build:** don't compile, build or run tests yourself unless the user
   explicitly allows it in this conversation; by default that is the
   user's job (`~/.claude/CLAUDE.md`). Never try to install a compiler.
   Say plainly the build is unrun and pending, under both infra layers,
   and never call a conversion verified on inspection alone.

### 10. Commit and push to `claude_<function-name>_containers` *(gated by the global git-commit preference, overridable per run — see below)*
   **Decide whether to run this step:**
   1. `--disable_git_commit` → skip; report modified files for manual
      commit. This run's override wins.
   2. `--enable_git_commit` → run this step. This run's override wins.
   3. Otherwise, read `~/.claude/preferences.json`'s `git_commit_and_push`:
      `"auto"` → run it; `"manual"`, missing file/key, or unparseable
      JSON → skip and report modified files. Unreadable configuration
      always means "manual"; never push on an unconfigured machine.

   Create (or check out) `claude_<lowercased_$1>_containers` from
   whatever branch the checkout is currently on. Stage every modified
   file and commit with a message naming the converted subroutine, the
   arguments that became containers, and the call sites updated. Use a
   HEREDOC so the trailer is preserved verbatim:

   ```
   BRANCH="claude_$(echo "$1" | tr '[:upper:]' '[:lower:]')_containers"
   git -C "$0" checkout -B "$BRANCH"
   git -C "$0" add -A
   git -C "$0" commit -m "$(cat <<'EOF'
   <one-line summary of the container conversion for $1>

   <1-3 line body: which dummies became containers, which grid-derived
   arrays were added, which call sites were updated, build status>

   Co-authored-by: Claude <noreply@anthropic.com>
   EOF
   )"
   git -C "$0" push -u origin "$BRANCH"
   ```

   If the push is rejected because the branch already exists upstream
   with unrelated history, stop and surface the conflict to the user
   rather than force-pushing.

## Versioning marker

Every Fortran file this skill creates or modifies gets a `!!SKILLS: 0.3.2`
marker line — the shared version number for this whole skill family,
not just this one skill (bump every skill file's `version:` field and
this marker in lockstep when any of them changes in a way that affects
generated code, never per-skill). If the file doesn't already have one,
add it as its own line immediately after the file's existing license/
header comment block, before the `module` statement. If it already has
one, update it to the current version rather than adding a second line.
Deliberately grep-able (`grep -rn "!!SKILLS:"`) and meant to be stripped
later, once these markers are no longer useful.

## Hard rules

- Never skip the `!!SKILLS: 0.3.2` marker on a file this skill touches,
  and never add a second marker line if one already exists — update it
  in place instead.
- Do not rename the subroutine, and do not create a `*_fortran` variant
  or a wrapper pair. That is `generate_cpp_bridge`'s job.
- Do not add a dispatcher, `select case (mode)`, `getenv_mode`, capture
  mode, `io_recorder`, `#ifdef _TIM`, `%to_c`, or a `bind(C)` interface
  where none exists — always `generate_cpp_bridge`'s job. If a target
  already has this content, convert normally (item 5) but never edit
  inside `select case…end select` except renaming a `CS%x`/`GV%x`
  reference, and never touch a `bind(C)` interface (lessons §10).
- Do not drop a `!<` doc comment — every dummy that had one keeps it
  word for word, including when a two-line declaration collapses to one
  (Step 4, where they're usually lost). Every added dummy gets a new one.
- Do not leave a declaration or the `subroutine <name>(...)` header over
  100 characters. Never fix it by shortening text or dropping an
  argument — wrap the comment or move the `&` continuation instead
  (Step 4 items 5–6, Step 9).
- Do not edit inside a loop body — only declarations and `%view` calls
  change; the math is untouched.
- Do not change the numerical result — a conversion is inert by
  construction; any difference is a bug.
- Do not pass a written container as `intent(out)` — always
  `intent(inout)` (lessons §9 #6).
- Do not skip `source=` for an `intent(out)` array as "pure output" —
  that reads garbage unless the write covers the array in full, which a
  `Box_t`-scoped write almost never does (lessons §6, §9 #2). Copy in by
  default; skip only when confirmed safe.
- Do not leave a `grow`/`shrink` result unfreed (lessons §9 #3).
- Do not use `%dup` — folded into `%alloc` (lessons §4.6).
- Do not convert a derived type's members (`BT_cont_type`) unless asked.
- Do not decide unilaterally whether an `optional` array dummy
  containerizes — Step 2 stops and asks per dummy; both outcomes are
  common.
- Do not leave a call site unupdated.
- Do not introduce a guarded `%view` (Step 5) without a `nullify` for
  that pointer in the same edit, even if every read sits inside the same
  guard (lessons §9 #16).
- Do not leave a container's `%associated()` as, or inside, an `!$omp
  target`/`teams`/`parallel`/`loop`-family construct's condition —
  precompute a scalar outside it (lessons §9 #17); a distinct bug from
  the `nullify` rule above, check both.
- Do not pass a container actual to a dummy Step 2 classified as raw, or
  vice versa (Step 8, Step 9).
- Do not attempt to install a Fortran compiler, or claim an unrun build
  passed.
- Do not add comments about the conversion itself — that a subroutine
  was migrated, why a dummy is now a pointer, that a `nullify` is
  mandatory (lessons §7). A new dummy's `!<` comment describes its
  physical meaning only — not an exception to this rule.
- Do not size anything off a `%view`-derived pointer's `LBOUND`/`UBOUND`
  (Step 2a) — only another dummy's own `%lb`/`%ub` is legal there; a
  `%view` pointer is only associated partway through the body, so
  referencing it in a specification expression is undefined behavior.
- Do not use an `optional` dummy's `%lb`/`%ub` to size another local or
  dummy (Step 2a) — an absent optional's container is disassociated, so
  its bounds would be too.
- Do not create a per-field container dummy, or a per-call
  `alloc(source=G%<f>)`, for a `G%` array. It comes from `Gcore`/`Gobc`
  (Step 2b). A field in neither structure is a question for the user.
- Do not create a new container dummy for a `G%`/`GV%` field when a
  local container already exists from that same source (Step 2b) —
  promote the existing local instead.
- Do not touch the body during a Step 2b promotion beyond the literal
  `G%<field>` → `<field>` substitution — the one narrow exception to the
  no-body-edits rule.
- Do not cherry-pick which qualifying fields to promote in a Step 2b
  pass — extract every distinct one the grep finds and confirm the full
  list with the user.

Anything not covered here: consult lessons §9 before improvising.

## Output to the user on success

Report:

1. The subroutine converted, and its new argument list.
2. Which dummies became containers, and which grid-derived arrays were
   added as new container dummies.
3. Which of `G` / `GV` / `US` / `CS` were dropped from the signature, if
   any.
4. Every call site updated (file and line).
5. Any callee still taking raw arrays, as a follow-up conversion
   candidate (Step 7).
6. **Hoisting candidates:** every container with a loop-invariant grid-
   metadata source, and the call sites that now rebuild it — worth a
   `hoist_container_marshalling` pass once this subroutine's
   conversions are complete and every call site is visible at once.
7. Every `optional` array dummy found, the decision for each (converted
   vs. left raw), and whether any converted one needed conditional
   allocation at a call site.
8. Build status — infra layers built, or explicitly why none were.
9. Committed, or the list of modified files for manual commit.
10. **Step 2a rewrites:** every local/left-raw dummy whose sizing was
    rewritten, which macro and existing dummy/dimension it now uses —
    plus every one considered but not rewired.
11. **Step 2b promotions:** every distinct field found, new dummy vs.
    merged-into-existing-local (name it), and the re-examined
    `G`/`GV`/`US`/`CS` drop decision. A caller needing its own Step 2b
    pass to supply a newly threaded dummy is a follow-up, not converted
    here.
