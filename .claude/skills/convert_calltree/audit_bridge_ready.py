#!/usr/bin/env python3
"""Bridge-readiness audit for a convert_calltree campaign (Phase 2 exit test, SKILL.md Step 6b).

Usage: audit_bridge_ready.py <file.F90> --targets a,b,... [--primitives c,d,...] [--forbid e,f,...]

  --targets     the Phase 3 bridge targets in the file
  --primitives  the per-point kernels the targets call (Step 2)
  --forbid      extra calls that must not appear in a target's body, e.g. calls into other
                modules that the plan routes through Fortran helpers (horizontal_viscosity:
                calc_QG_slopes,calc_QG_Leith_viscosity,thickness_to_dz,barotropic_get_tav,
                thickness_diffuse_get_KH,ZB2020_copy_gradient_and_thickness,ZB2020_lateral_stress)

Example (vertvisc_family):
  audit_bridge_ready.py src/parameterizations/vertical/MOM_vert_friction.F90 \
      --targets vertvisc_TR,vertvisc_coef_TR,vertvisc_remnant_TR,vertvisc_limit_vel \
      --primitives find_coupling_coef_k,find_coupling_coef,find_coupling_coef_gl90

Checks, per Phase 3 target:
  1. every dummy can cross bind(C) (container/box, non-pointer scalar, grid_core_type,
     a derived type with a bind(C) mirror AND a converter, the OBC pointer, or one opaque
     host-context `type(c_ptr) :: host` dummy that reaches Fortran-only helpers);
  2. no Fortran-only work in the body (diagnostics, halo passes incl. pass_var/pass_vector,
     checksums, reporting, any --forbid call, %alloc/%free/%rebound of a dummy or a component of one
     (local scratch containers are allowed),
     MOM_error, I/O, G%/GV%/US% references, SZ*_ macros, container %alloc/%free);
  3. every do concurrent iterates over a Box_t;
and, for the per-column primitives:
  4. no tv/VarMix/G/GV/US dummies (OBC pointer allowed).
Prints PASS/FAIL per item and exits non-zero if anything fails.
"""
import re
import sys
import argparse

ap = argparse.ArgumentParser()
ap.add_argument('file')
ap.add_argument('--targets', required=True)
ap.add_argument('--primitives', default='')
ap.add_argument('--forbid', default='')
opts = ap.parse_args()
src = open(opts.file).read()
TARGETS = [t for t in opts.targets.split(',') if t]
PRIMITIVES = [p for p in opts.primitives.split(',') if p]

mirrors = {m.lower() for m in re.findall(r'^type, bind\(C\) :: (\w+)', src, re.M)}
converters = {}
for m in re.finditer(r'^function (\w+_to_c)\((\w+)\) result\(\w+\)\n((?:.*\n){1,3})', src, re.M):
    t = re.search(r'type\((\w+)\)', m.group(3))
    if t: converters[t.group(1).lower()] = m.group(1)


def sub(name):
    m = re.search(r'^(?:pure )?subroutine ' + name + r'\(.*?^end subroutine ' + name + r'\s*$',
                  src, re.M | re.S)
    return m.group(0) if m else None


def stmts(t):
    L = t.split('\n'); out = []; i = 0
    while i < len(L):
        st = L[i] if L[i].lstrip().startswith('!$') else re.sub(r'!.*', '', L[i]); j = i
        while st.rstrip().endswith('&'):
            j += 1; st = st.rstrip()[:-1] + ' ' + re.sub(r'!.*', '', L[j]).strip().lstrip('&')
        out.append(st); i = j + 1
    return out


def dummy_types(t):
    st = stmts(t)
    names = [a.strip() for a in re.search(r'\((.*?)\)', st[0]).group(1).split(',') if a.strip()]
    decl = {}
    for x in st[1:]:
        if '::' in x and not x.lstrip().startswith('!'):
            typ = re.sub(r'\s+', ' ', x.split('::')[0].strip())
            b = x.split('::', 1)[1]
            while re.search(r'\([^()]*\)', b): b = re.sub(r'\([^()]*\)', '', b)
            for v in b.split(','):
                decl.setdefault(v.split('=')[0].strip().lower(), typ)
    return [(n, decl.get(n.lower(), '?')) for n in names]


def mirror_of(tn):
    tn = tn.lower()
    cands = {tn + '_c', tn.replace('_type', '_c'), tn.replace('_cs', '_c')}
    return next((c for c in cands if c in mirrors), None)


def classify(name, typ):
    tl = typ.lower()
    m = re.match(r'type\((\w+)\)', tl)
    if 'optional' in tl and not (m and m.group(1) == 'ocean_obc_type'):
        return False, 'optional dummy'
    if m:
        tn = m.group(1)
        if tn in ('realarray_t', 'intarray_t', 'logicalarray_t', 'box_t'): return True, 'container/box'
        if tn == 'grid_core_type': return True, 'grid_core_C'
        if tn == 'ocean_obc_type' and 'pointer' in tl: return True, 'OBC -> c_ptr'
        if tn == 'c_ptr' and name.lower() == 'host': return True, 'host context -> c_ptr'
        mc = mirror_of(tn)
        if mc and tn in converters: return True, 'mirror %s + %s' % (mc, converters[tn])
        if mc: return False, 'mirror %s but no converter' % mc
        return False, 'no bind(C) mirror for %s' % tn
    if re.match(r'(real|integer|logical)\b', tl):
        if 'dimension' in tl or 'pointer' in tl: return False, 'raw array/pointer'
        return True, 'scalar'
    if tl.startswith('character'): return False, 'character'
    return False, 'unrecognized type'


EXTRA = [re.escape(x.strip()) for x in opts.forbid.split(',') if x.strip()]
FORBIDDEN = (r'\b(post_data|post_product\w*|query_averaging_enabled|create_group_pass|do_group_pass|'
             r'pass_var|pass_vector|'
             r'\w*chksum|write_[uv]_accel|MOM_error|MOM_read_data|open|write' +
             ''.join('|' + x for x in EXTRA) + r')\s*\(')

ok_all = True
print('%-22s %-6s %-6s %-6s  details' % ('target', 'args', 'body', 'loops'))
for t in TARGETS:
    s = sub(t)
    if s is None:
        print('%-22s MISSING' % t); ok_all = False; continue
    bad_args = [(n, why) for n, typ in dummy_types(s) for ok, why in [classify(n, typ)] if not ok]
    body = '\n'.join(stmts(s)[1:])
    code = '\n'.join(l for l in body.split('\n') if not l.lstrip().startswith('!$'))
    bad_body = sorted(set(re.findall(FORBIDDEN, code)))
    bad_body += sorted(set('%s%%' % x for x in re.findall(r'(?<![%\w])(G|GV|US)%\w+', code)))
    bad_body += sorted(set(re.findall(r'SZ[IJK]B?_\(\w+\)', code)))
    # %alloc/%free/%rebound are Fortran-only, but only matter on data that crosses the bridge or
    # persists: a dummy or a component of one.  Local scratch inside the Fortran implementation is
    # allowed (the C++ side manages its own; CorAdCalc_TR precedent).
    dummy_names = {n.lower() for n, _ in dummy_types(s)}
    for base, op in re.findall(r'(?<![%\w])(\w+)(?:%\w+)*%(alloc|free|rebound)\s*\(', code):
        if base.lower() in dummy_names:
            bad_body.append('%s%%...%%%s()' % (base, op))
    bad_body = sorted(set(bad_body))
    heads = re.findall(r'do concurrent \(([^\n]*)', code)
    # A loop clipped by OBC segment bounds runs only with OBC active, which the C++ aborts on.
    bad_loops = [h[:60] for h in heads if 'OBC%' not in h
                 if not all(re.match(r'\s*\w+\s*=\s*\w+%idxS\(\d\):\w+%idxE\(\d\)\s*$', r)
                            for r in re.split(r',(?![^()]*\))', h.split(')')[0] if False else
                                              re.sub(r'\)\s*(DO_LOCALITY.*)?$', '', h))
                            if '=' in r and ':' in r)]
    a, b, c = not bad_args, not bad_body, not bad_loops
    ok_all &= a and b and c
    print('%-22s %-6s %-6s %-6s' % (t, 'PASS' if a else 'FAIL', 'PASS' if b else 'FAIL',
                                    'PASS' if c else 'FAIL'))
    for n, why in bad_args: print('    arg  %-15s %s' % (n, why))
    for x in bad_body: print('    body %s' % x)
    for x in bad_loops: print('    loop %s' % x)

print()
for p in PRIMITIVES:
    s = sub(p)
    if s is None:
        print('%-24s MISSING' % p); ok_all = False; continue
    bad = [(n, typ) for n, typ in dummy_types(s)
           if re.search(r'type\((thermo_var_ptrs|varmix_cs|ocean_grid_type|verticalgrid_type|'
                        r'unit_scale_type)\)', typ.lower())]
    ok_all &= not bad
    print('%-24s %s %s' % (p, 'PASS' if not bad else 'FAIL', ' '.join('%s:%s' % x for x in bad)))

print('\nBRIDGE-READY' if ok_all else '\nNOT BRIDGE-READY')
sys.exit(0 if ok_all else 1)
