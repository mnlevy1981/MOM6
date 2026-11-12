#!/usr/bin/env python3
"""Classify the serial `do` loop nests of MOM6 subroutines for convert_serial_loops_to_do_concurrent.

Usage: check_loop_independence.py <file.F90> <subroutine>[,<subroutine>...]

For every serial nest (`do a=... ; do b=...` chained on one line, or nested on following lines
with nothing between the `do` statements) prints one of:
  MAP        every iteration writes only its own point: convertible to one `do concurrent`;
             lists the scalar temporaries that must go in DO_LOCALITY(local(...))
  REDUCTION  per-point accumulation over the outermost index: convertible by interchange
             (`do concurrent` over the inner indices, the reduction index serial inside)
  SEQUENTIAL with the reasons; never converted as a whole
Every serial `do` left outside a MAP or REDUCTION nest (the levels of a SEQUENTIAL nest, perfect or
not) then gets a verdict of its own:
  LEVEL-CONCURRENT  iterations over this one index are independent, its body (inner loops
             included) unchanged: it can become `do concurrent` where it stands; lists its
             `local` scalars and the passing level directly inside it that it merges with
  LEVEL-SERIAL      with the reasons; stays a serial `do`
It is deliberately conservative: anything it cannot prove is SEQUENTIAL.  Assumes distinct
array names never alias (true for distinct `%view` pointers of distinct containers and for
distinct dummies).  Derived-type designators (`CS%Kh_bg_xy(I,J)`, `G%IdxCu(I,j)`) are treated as
arrays; a write to a scalar component inside a nest is shared by all iterations, so SEQUENTIAL.
A nest with an `!$omp`/`!$acc` directive on it or inside it is SEQUENTIAL until the plan's
device-directive policy says otherwise."""
import re
import sys

KEYWORDS = {'if', 'then', 'else', 'elseif', 'endif', 'end', 'do', 'enddo', 'cycle', 'and', 'or',
            'not', 'true', 'false', 'abs', 'sqrt', 'max', 'min', 'merge', 'sign', 'atan', 'exp',
            'log', 'real', 'int', 'nint', 'mod', 'huge', 'tiny', 'epsilon', 'present', 'associated'}
FORBIDDEN = r'\b(exit|return|goto|stop|print|write|read|open|close)\b'


def code(l):
    return re.sub(r'!.*', '', l)


def statements(lines):
    """(first line index, joined text) per statement, splitting on ';' and joining '&'."""
    out, cur, start = [], '', None
    for n, l in enumerate(lines):
        c = code(l).strip()
        if not c:
            continue
        if start is None:
            start = n
        if cur:
            c = c.lstrip('&').strip()
        if c.endswith('&'):
            cur += c[:-1].rstrip() + ' '
            continue
        for part in (cur + c).split(';'):
            if part.strip():
                out.append((start, part.strip()))
        cur, start = '', None
    return out


def refs(expr):
    """Array references name(args) at depth 0 of their own parentheses, plus bare names."""
    arr, i = [], 0
    for m in re.finditer(r'(?<![%\w])((?:[A-Za-z]\w*%)*[A-Za-z]\w*)\s*\(', expr):
        name, p, d = m.group(1), m.end(), 1
        q = p
        while q < len(expr) and d:
            d += {'(': 1, ')': -1}.get(expr[q], 0)
            q += 1
        arr.append((name, [a.strip() for a in split_args(expr[p:q - 1])]))
    bare = set(re.findall(r'(?<![%\w])([A-Za-z]\w*)\b(?!\s*\()', expr))
    return arr, bare


def split_args(s):
    a, d, c = [], 0, ''
    for ch in s:
        if ch == ',' and d == 0:
            a.append(c)
            c = ''
            continue
        d += {'(': 1, ')': -1}.get(ch, 0)
        c += ch
    a.append(c)
    return a


def analyse(nest_vars, body, declared_arrays):
    lv = [v.lower() for v in nest_vars]
    alias = {}        # scalar -> loop var it is an affine shift of (k = kstart + kk - 1)
    scalars_set, scalars_read_first, reasons = [], [], []
    writes, reads = {}, {}
    for _, st in body:
        low = st.lower()
        if re.search(FORBIDDEN, low):
            reasons.append('forbidden statement: ' + st[:50])
        if re.match(r'call\b', low):
            reasons.append('call (only pure/elemental allowed; check by hand): ' + st[:50])
        if re.match(r'(do|enddo|end do)\b', low):
            reasons.append('inner loop: ' + st[:50])
            continue
        m = re.match(r'(?:if\s*\((.*)\)\s*)?((?:[A-Za-z]\w*%)*[A-Za-z]\w*)(\s*\(([^=]*)\))?\s*=(?!=)\s*(.*)$', st)
        cond_and_rhs = st
        if m and not re.match(r'(if|else|elseif|end|endif|cycle)\b', m.group(2).lower()):
            tgt, idx, rhs = m.group(2), m.group(4), m.group(5)
            cond_and_rhs = (m.group(1) or '') + ' ' + rhs + ' ' + (idx or '')
            if idx is None and '%' in tgt:
                reasons.append('writes a derived-type component shared by all iterations: ' + tgt)
            elif idx is None:
                t = tgt.lower()
                if t in lv:
                    reasons.append('assigns a loop index')
                if t not in scalars_set:
                    scalars_set.append(t)
                for v in lv:
                    if re.fullmatch(r'\s*\w+\s*\+\s*%s\s*-\s*1\s*' % v, rhs, re.I) or \
                       re.fullmatch(r'\s*%s\s*[+-]\s*\w+\s*' % v, rhs, re.I):
                        alias[t] = v
            else:
                writes.setdefault(tgt.lower(), set()).add(tuple(a.strip().lower() for a in split_args(idx)))
        arr, bare = refs(cond_and_rhs)
        for name, args in arr:
            if name.lower() in declared_arrays or '%' in name:
                reads.setdefault(name.lower(), set()).add(tuple(a.lower() for a in args))
        for b in bare:
            b = b.lower()
            if b in KEYWORDS or b in lv:
                continue
            if b not in scalars_set and b not in scalars_read_first:
                scalars_read_first.append(b)
    carried = [s for s in scalars_set if s in scalars_read_first and s not in alias]
    if carried:
        reasons.append('scalar read before assignment (carried): ' + ', '.join(carried))

    def canon(t):
        return tuple(alias.get(x, x) for x in t)
    reduction = []
    for a, idxs in writes.items():
        for t in idxs:
            ct = canon(t)
            used = {x for x in ct if x in lv}
            if not all(x in lv or x == ':' for x in ct):
                reasons.append('write %s%s not indexed by the loop indices' % (a, t))
            elif used != set(lv):
                missing = [v for v in lv if v not in used]
                if missing == [lv[0]] and len(idxs) == 1 and canon(next(iter(reads.get(a, {t})))) == ct:
                    reduction.append(a)
                else:
                    reasons.append('write %s%s misses loop index %s' % (a, t, ','.join(missing)))
            for r in reads.get(a, ()):
                if canon(r) != ct:
                    reasons.append('%s written at %s but read at %s' % (a, t, r))
    if reasons:
        return 'SEQUENTIAL', reasons
    local = [s for s in scalars_set if s not in lv]
    if reduction:
        return 'REDUCTION', ['accumulates %s over %s' % (', '.join(reduction), lv[0]),
                             'local: ' + (', '.join(local) or '-')]
    return 'MAP', ['local: ' + (', '.join(local) or '-')]



def header(st):
    """('serial', [var], bounds, step) or ('concurrent', [vars], bounds, locals) for a `do`, else None."""
    m = re.match(r'do\s+concurrent\s*\(', st, re.I)
    if m:
        p, d = m.end(), 1
        q = p
        while q < len(st) and d:
            d += {'(': 1, ')': -1}.get(st[q], 0)
            q += 1
        specs = split_args(st[p:q - 1])
        names = [x.split('=')[0].strip().lower() for x in specs]
        bounds = ' '.join(x.split('=', 1)[1] for x in specs if '=' in x)
        return 'concurrent', names, bounds, st[q:]
    m = re.match(r'do\s+(\w+)\s*=\s*(.*)$', st, re.I)
    if m:
        parts = split_args(m.group(2))
        return 'serial', [m.group(1).lower()], ' '.join(parts[:2]), (parts[2].strip() if len(parts) > 2 else None)
    return None


def loop_end(sts, t):
    """Index of the `enddo` closing the loop opened at statement t."""
    depth, b = 0, t
    while b < len(sts):
        low = sts[b][1].lower()
        if re.match(r'do\b', low):
            depth += 1
        elif re.match(r'(enddo|end\s*do)\b', low):
            depth -= 1
            if depth == 0:
                return b
        b += 1
    raise ValueError('unclosed loop at statement %d' % t)


def level(var, step, body):
    """Can the serial loop over `var` become `do concurrent` with `body` (inner loops included)
    unchanged?  Every array the body writes must take `var` as a bare subscript at one fixed
    position in every reference to it, so distinct iterations touch disjoint elements."""
    v = var.lower()
    reasons = []
    if step is not None:
        reasons.append('loop step ' + step)
    scalars_set, read_first, scoped, inner = [], [], set(), []
    written, allrefs, bare_names = set(), [], set()
    for n, (_, st) in enumerate(body):
        low = st.lower()
        if re.search(FORBIDDEN, low):
            reasons.append('forbidden statement: ' + st[:50])
        if re.match(r'call\b', low):
            reasons.append('call (only pure/elemental allowed; check by hand): ' + st[:50])
        if re.match(r'(enddo|end\s*do)\b', low):
            continue
        h = header(st)
        text = st
        if h:
            kind, names, bounds, _ = h
            text = bounds
            if kind == 'concurrent':
                scoped.update(names)      # construct entities of the inner do concurrent
            else:
                for x in names:
                    if x == v:
                        reasons.append('assigns the loop index')
                    else:
                        inner.append((n, x))   # local, unless that inner loop is converted too
        else:
            m = re.match(r'(?:if\s*\((.*)\)\s*)?((?:[A-Za-z]\w*%)*[A-Za-z]\w*)(\s*\(([^=]*)\))?\s*=(?!=)\s*(.*)$', st)
            if m and not re.match(r'(if|else|elseif|end|endif|cycle)\b', m.group(2).lower()):
                tgt, idx, rhs = m.group(2), m.group(4), m.group(5)
                text = (m.group(1) or '') + ' ' + rhs + ' ' + (idx or '')
                if idx is None and '%' in tgt:
                    reasons.append('writes a derived-type component shared by all iterations: ' + tgt)
                elif idx is None:
                    t = tgt.lower()
                    if t == v:
                        reasons.append('assigns the loop index')
                    # the right-hand side is read before the assignment happens
                    arr, bare = refs(text)
                    for name, args in arr:
                        allrefs.append((name.lower(), tuple(a.lower() for a in args)))
                    for b in bare:
                        b = b.lower()
                        if b not in KEYWORDS and b != v and b not in scoped and b not in scalars_set \
                                and b not in read_first and b not in [x for _, x in inner]:
                            read_first.append(b)
                    bare_names.update(x.lower() for x in bare)
                    if t not in scalars_set:
                        scalars_set.append(t)
                    continue
                else:
                    written.add(tgt.lower())
                    allrefs.append((tgt.lower(), tuple(a.strip().lower() for a in split_args(idx))))
        arr, bare = refs(text)
        for name, args in arr:
            allrefs.append((name.lower(), tuple(a.lower() for a in args)))
        bare_names.update(x.lower() for x in bare)
        for b in bare:
            b = b.lower()
            if b in KEYWORDS or b == v or b in scoped or b in [x for _, x in inner]:
                continue
            if b not in scalars_set and b not in read_first:
                read_first.append(b)
    carried = [x for x in dict.fromkeys(scalars_set + [x for _, x in inner]) if x in read_first]
    if carried:
        reasons.append('scalar read before assignment (carried): ' + ', '.join(carried))
    pos = {}
    for a, t in allrefs:
        if a not in written:
            continue
        at = [p for p, x in enumerate(t) if x == v]
        other = [p for p, x in enumerate(t) if x != v and re.search(r'\b%s\b' % v, x)]
        if len(at) != 1 or other:
            reasons.append('%s%s does not take %s as one bare subscript' % (a, t, v))
        elif pos.setdefault(a, at[0]) != at[0]:
            reasons.append('%s%s takes %s at another position' % (a, t, v))
    for a in written & bare_names:
        reasons.append('whole-array reference to written %s' % a)
    return list(dict.fromkeys(reasons)), [x for x in scalars_set if x != v and x not in scoped], inner


def h_var(st):
    return header(st)[1][0]


def main():
    path, names = sys.argv[1], [n.lower() for n in sys.argv[2].split(',')]
    L = open(path).read().split('\n')
    i = 0
    while i < len(L):
        m = re.match(r'\s*(?:pure\s+|elemental\s+)*subroutine\s+(\w+)', code(L[i]), re.I)
        if not m or m.group(1).lower() not in names:
            i += 1
            continue
        s = i
        e = next(j for j in range(s, len(L)) if re.match(r'\s*end\s+subroutine\b', code(L[j]), re.I))
        sub = L[s:e + 1]
        declared = set()
        for _, st in statements(sub):
            if '::' in st and ('dimension' in st.lower() or 'pointer' in st.lower() or '(' in st.split('::')[1]):
                rhs = st.split('::', 1)[1]
                if 'dimension' in st.lower().split('::')[0] or 'pointer' in st.lower().split('::')[0]:
                    names_ = [n.split('(')[0].split('=')[0].strip().lower() for n in split_args(rhs)]
                else:
                    names_ = [n.split('(')[0].strip().lower() for n in split_args(rhs) if '(' in n]
                declared.update(names_)
        sts = statements(sub)
        covered, serial_idx = set(), set()
        k = 0
        while k < len(sts):
            ln, st = sts[k]
            if not re.match(r'do\s+\w+\s*=', st, re.I) or re.match(r'do\s+concurrent', st, re.I):
                k += 1
                continue
            nest, j = [], k
            while j < len(sts) and re.match(r'do\s+(\w+)\s*=', sts[j][1], re.I) and \
                    not re.match(r'do\s+concurrent', sts[j][1], re.I):
                nest.append(re.match(r'do\s+(\w+)', sts[j][1], re.I).group(1))
                j += 1
            depth, b = len(nest), j
            while depth:
                t = sts[b][1].lower()
                if re.match(r'do\b', t):
                    depth += 1
                elif re.match(r'(enddo|end do)\b', t):
                    depth -= 1
                b += 1
            body = sts[j:b - len(nest)] if all(re.match(r'(enddo|end do)', sts[x][1], re.I)
                                               for x in range(b - len(nest), b)) else sts[j:b - 1]
            kind, why = analyse(nest, body, declared)
            last = sts[b - 1][0]
            before = ln - 1
            while before >= 0 and not sub[before].strip():
                before -= 1
            span = sub[max(before, 0):last + 1]
            if any(re.match(r'\s*!\$(omp|acc)\b', x, re.I) for x in span):
                kind, why = 'SEQUENTIAL', ['OpenMP/OpenACC directive on or in the nest: apply the '
                                           "plan's device-directive policy"] + (why if kind == 'SEQUENTIAL' else [])
            print('%-26s line %5d  %-11s do %-18s %s' % (m.group(1), s + ln + 1, kind, ','.join(nest),
                                                        '; '.join(why)))
            if kind != 'SEQUENTIAL':
                covered.update(range(k, b))
            if kind == 'REDUCTION':
                serial_idx.add(k)          # the reduction index stays a serial do, moved inside
            k = j if kind == 'SEQUENTIAL' else b   # inner nests of a sequential loop are analysed too
        verdict = {}
        for t, (ln, st) in enumerate(sts):
            h = header(st)
            if t in covered or not h or h[0] != 'serial':
                continue
            end = loop_end(sts, t)
            why, local, inner = level(h[1][0], h[3], sts[t + 1:end])
            before = ln - 1
            while before >= 0 and not sub[before].strip():
                before -= 1
            if any(re.match(r'\s*!\$(omp|acc)\b', x, re.I) for x in sub[max(before, 0):sts[end][0] + 1]):
                why = ["OpenMP/OpenACC directive on or in the loop: apply the plan's device-directive policy"] + why
            verdict[t] = (end, not why, why, local, [(t + 1 + n, x) for n, x in inner])
        for t, (end, ok, why, local, inner) in sorted(verdict.items()):
            ln = sts[t][0]
            if ok:
                # an inner loop that stays serial (or is a REDUCTION's serial index) defines its index
                # in every iteration: local; one converted to do concurrent scopes its own index
                stays = [x for u, x in inner if (u in verdict and not verdict[u][1]) or u in serial_idx]
                local = list(dict.fromkeys(local + stays))
                inn = verdict.get(t + 1)
                if inn and inn[1] and inn[0] == end - 1:
                    local = [x for x in local if x != h_var(sts[t + 1][1])]
                merge = ''
                if inn and inn[1] and inn[0] == end - 1:
                    merge = '; merges with the level at line %d' % (s + sts[t + 1][0] + 1)
                print('%-26s line %5d  LEVEL-CONCURRENT do %-9s local: %s%s' % (
                    m.group(1), s + ln + 1, h_var(sts[t][1]), ', '.join(local) or '-', merge))
            else:
                print('%-26s line %5d  LEVEL-SERIAL     do %-9s %s' % (m.group(1), s + ln + 1, h_var(sts[t][1]),
                                                                       '; '.join(why)))
        i = e + 1


if __name__ == '__main__':
    main()
