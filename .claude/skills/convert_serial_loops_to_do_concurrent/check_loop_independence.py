#!/usr/bin/env python3
"""Classify the serial `do` loop nests of MOM6 subroutines for convert_serial_loops_to_do_concurrent.

Usage: check_loop_independence.py <file.F90> <subroutine>[,<subroutine>...]

For every serial nest (`do a=... ; do b=...` chained on one line, or nested on following lines
with nothing between the `do` statements) prints one of:
  MAP        every iteration writes only its own point: convertible to one `do concurrent`;
             lists the scalar temporaries that must go in DO_LOCALITY(local(...))
  REDUCTION  per-point accumulation over the outermost index: convertible by interchange
             (`do concurrent` over the inner indices, the reduction index serial inside)
  SEQUENTIAL with the reasons; never converted
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
            k = j if kind == 'SEQUENTIAL' else b   # inner nests of a sequential loop are analysed too
        i = e + 1


if __name__ == '__main__':
    main()
