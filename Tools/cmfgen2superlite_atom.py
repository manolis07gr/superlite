#!/usr/bin/env python3
# This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
"""
Convert a CMFGEN model atom (Hillier & Miller 1998; atomic data distributed
with CMFGEN at https://sites.pitt.edu/~hillier/web/CMFGEN.htm) into a
SuperLite NLTE model-atom file.

Input (one ion):
  --osc   oscillator file (energy levels, f-values, A-values)
  --phot  photoionization file (route A: to the ground state of the next ion)
  --phot-b  (repeatable) files of the further routes (B, C, ...: to excited
          terms of the next ion), with --next-osc, the oscillator file of the
          next ion (energy of the final term: g-weighted mean of its levels,
          matched by name as in CMFGEN's RDPHOT_GEN)
  --col   collision-strength file (Omega(T), T in 10^4 K)
  --hyd   directory with hyd_l_data.dat and gbf_n_data.dat (hydrogenic
          cross-sections; needed for cross-section types 2, 3 and 8)
  --nlev  number of levels to keep (lowest in the osc file)
  --ftos  optional super-level file; levels then carry their super-level
          index (only for information: SuperLite treats every level)

The photoionization cross-sections are evaluated with the formulae of
CMFGEN's SUB_PHOT_GEN (types 1-9, 20, 21) and tabulated against
u = nu/nu_edge.  Further routes are added to route A (one cross-section to
the ground state of the next ion): with the final term in Boltzmann
equilibrium with the ground state of the next ion at T_e (metastable terms
of dense gas), the Milne recombination of the sum is that of the separate
routes.  Collision strengths are copied; SuperLite computes the
missing ones at run time with CMFGEN's rules (OMEGA_GEN_V3): van Regemorter
with the Auer & Mihalas (neutrals) or max(0.2|0.7, 0.276 e^x E1(x)) (ions)
Gaunt factor for f > 1e-5, OMEGA_SET otherwise; Seaton's formula for
collisional ionization.  The flag same_n of each line (CMFGEN SAME_N) is
written for the 0.7 Gaunt factor of Delta n = 0 transitions of ions.

Run with python3 -I on downloaded data.
"""
import argparse
import datetime
import math
import os
import re
import sys

import numpy as np

CLIGHT = 2.99792458e10
RINF = 109737.31          # cm^-1
MBARN = 1e-18
EV_TO_HZ15 = 0.241798840766


def read_lines(path):
    with open(path, 'r', errors='replace') as f:
        return f.read().splitlines()


def header_value(lines, key, conv=float):
    for s in lines:
        if key in s:
            return conv(s.split('!')[0].split()[0].replace('D', 'E').replace('d', 'e'))
    raise ValueError('key not found: ' + key)


def fnum(s):
    return float(s.replace('D', 'E').replace('d', 'e'))


# ---------------------------------------------------------------------------
# oscillator file
# ---------------------------------------------------------------------------
def read_osc(path):
    lines = read_lines(path)
    # the header block may be repeated; take the first numeric header
    nlev = header_value(lines, '!Number of energy levels', lambda s: int(float(s)))
    chi = header_value(lines, '!Ionization energy')
    zion = header_value(lines, '!Screened nuclear charge')
    ntr = header_value(lines, '!Number of transitions', lambda s: int(float(s)))
    # levels start after the last '!Number of transitions' header line
    i0 = max(i for i, s in enumerate(lines) if '!Number of transitions' in s) + 1
    levels = []
    i = i0
    while len(levels) < nlev:
        s = lines[i]
        i += 1
        t = s.split()
        if len(t) < 3:
            continue
        try:
            g = fnum(t[1])
            e = fnum(t[2])
        except ValueError:
            continue
        # ID: integer column after lambda (6th entry) if present
        lid = int(float(t[6])) if len(t) > 6 else len(levels) + 1
        levels.append(dict(name=t[0], g=g, e=e, id=lid))
    # transitions: lines ending with "i- j  [n]" (2023 format: then |-separated
    # columns); some files list the transitions twice: the first entry is kept
    pat = re.compile(r'\s(\S+)\s+(\S+)\s+(-?\S+)\s+(\d+)\s*-\s*(\d+)(?:\s+(\d+))?\s*$')
    trans = []
    seen = set()
    for s in lines[i:]:
        s = s.split('|')[0].rstrip()          # 2023 format: observed lambda, accuracy
        m = pat.search(s)
        if not m:
            continue
        try:
            f = fnum(m.group(1))
            a = fnum(m.group(2))
            lam = fnum(m.group(3))
        except ValueError:
            continue
        lo, up = int(m.group(4)), int(m.group(5))
        if (lo, up) in seen:
            continue
        seen.add((lo, up))
        trans.append(dict(l=lo, u=up, f=f, a=a, lam=lam))
        if len(trans) == ntr:
            break                           # CMFGEN reads the stated number
    if len(trans) != ntr:
        print('warning: %s: %d transitions read, header says %d' % (path, len(trans), ntr),
              file=sys.stderr)
    return dict(nlev=nlev, chi=chi, zion=zion, levels=levels, trans=trans)


# ---------------------------------------------------------------------------
# super-level links (RD_F_TO_S_IDS_V3)
# ---------------------------------------------------------------------------
def read_ftos(path, nlev):
    lines = read_lines(path)
    entry = header_value(lines, '!Entry number of link', lambda s: int(float(s)))
    i0 = max(i for i, s in enumerate(lines) if '!Entry number of link' in s) + 1
    sl = []
    for s in lines[i0:]:
        t = s.split()
        if not t or t[0].startswith('!'):
            continue
        sl.append(int(float(t[entry - 1])))
        if len(sl) == nlev:
            break
    # compact numbering 1..ns in order of first appearance
    remap = {}
    out = []
    for k in sl:
        if k not in remap:
            remap[k] = len(remap) + 1
        out.append(remap[k])
    return out


# ---------------------------------------------------------------------------
# hydrogenic tables (RD_HYD_BF_DATA)
# ---------------------------------------------------------------------------
class Hyd:
    def __init__(self, d):
        self.l_ok = self.n_ok = False
        if d is None:
            return
        pl = os.path.join(d, 'hyd_l_data.dat')
        pn = os.path.join(d, 'gbf_n_data.dat')
        if os.path.exists(pl):
            self._read_l(pl)
        if os.path.exists(pn):
            self._read_n(pn)

    @staticmethod
    def _nums(lines):
        out = []
        for s in lines:
            for t in s.split():
                out.append(fnum(t))
        return out

    def _read_l(self, p):
        lines = read_lines(p)
        k = next(i for i, s in enumerate(lines) if 'Maximum principal quantum number' in s)
        self.lmaxn = int(fnum(lines[k].split()[0]))
        self.nperl = int(fnum(lines[k + 1].split()[0]))
        self.ldel = fnum(lines[k + 3].split()[0])
        vals = self._nums(lines[k + 4:])
        self.lcross = {}
        p = 0
        for n in range(1, self.lmaxn + 1):
            for l in range(n):
                rn, rl, cnt = int(vals[p]), int(vals[p + 1]), int(vals[p + 2])
                assert rn == n and rl == l and cnt == self.nperl, (n, l, rn, rl, cnt)
                p += 3
                self.lcross[(n, l)] = np.array(vals[p:p + cnt])
                p += cnt
        self.l_ok = True

    def _read_n(self, p):
        lines = read_lines(p)
        k = next(i for i, s in enumerate(lines) if 'Maximum principal quantum number' in s)
        self.nmaxn = int(fnum(lines[k].split()[0]))
        self.npern = int(fnum(lines[k + 1].split()[0]))
        self.ndel = fnum(lines[k + 3].split()[0])
        vals = self._nums(lines[k + 4:])
        self.ngaunt = {}
        p = 0
        for n in range(1, self.nmaxn + 1):
            rn, cnt = int(vals[p]), int(vals[p + 1])
            assert rn == n and cnt == self.npern
            p += 2
            self.ngaunt[n] = np.array(vals[p:p + cnt])
            p += cnt
        self.n_ok = True

    def log_cross_l(self, n, l, u):
        """log10 of the H (n,l) cross-section in program units (cgs*1e10)"""
        tab = self.lcross[(n, l)]
        rj = math.log10(u) / self.ldel
        j = int(rj)
        t1 = rj - j
        j += 1                       # 1-based
        if j < self.nperl:
            return t1 * tab[j] + (1 - t1) * tab[j - 1]
        return (tab[-1] - tab[-2]) * (rj - self.nperl) + tab[-1]

    def gaunt_n(self, n, u):
        if n > 30:
            return 1.0
        tab = self.ngaunt[n]
        rj = math.log10(u) / self.ndel
        j = int(rj)
        t1 = rj - j
        j += 1
        if j < self.npern:
            return t1 * tab[j] + (1 - t1) * tab[j - 1]
        t1 = math.log10(tab[-2] / tab[-1])
        return tab[-1] * 10 ** (t1 * (self.npern - rj))


# ---------------------------------------------------------------------------
# photoionization file (RDPHOT_GEN_V2) and cross-sections (SUB_PHOT_GEN)
# ---------------------------------------------------------------------------
def read_phot(path):
    lines = read_lines(path)
    hdr = {}
    for key in ['!Number of energy levels', '!Number of photoionization routes',
                '!Screened nuclear charge', '!Excitation energy of final state',
                '!Statistical weight of ion', '!Split J levels',
                '!Sigma used for Gaussian smoothing (km/s)']:
        for s in lines:
            if key in s:
                hdr[key] = s.split('!')[0].split()[0]
                break
    terms = {}
    i = 0
    n = len(lines)
    while i < n:
        s = lines[i]
        if '!Configuration name' in s:
            name = s.split('!')[0].strip()
            j = i + 1
            while '!Type of cross-section' not in lines[j]:
                j += 1
            typ = int(fnum(lines[j].split()[0]))
            npts = int(fnum(lines[j + 1].split()[0]))
            vals = []
            k = j + 2
            need = npts if typ < 20 else 2 * npts
            while len(vals) < need:
                for t in lines[k].split():
                    if t.startswith('!'):
                        break
                    vals.append(fnum(t))
                k += 1
            if typ < 20:
                terms[name] = (typ, np.array(vals[:npts]), None)
            else:
                a = np.array(vals[:2 * npts]).reshape(npts, 2)
                terms[name] = (typ, a[:, 1].copy(), a[:, 0].copy())
            i = k
            continue
        i += 1
    return hdr, terms


def norm_name(s):
    return s.split('[')[0].replace('_', '').replace(' ', '')


def final_term_energy(final, levels):
    """g-weighted energy [cm^-1] of the levels of the next ion whose name
    (without J) is one of the '/'-separated alternatives of 'final'."""
    alts = [norm_name(a) for a in final.split('/')]
    gs = es = 0.0
    for lv in levels:
        if norm_name(lv['name']) in alts:
            gs += lv['g']
            es += lv['g'] * lv['e']
    return es / gs if gs > 0 else None


def term_for_level(levname, terms, split_j):
    key = levname
    if not split_j:
        k = levname.find('[')
        if k > 0:
            key = levname[:k]
    if key in terms:
        return terms[key]
    # alternate names separated by '/'
    for t, v in terms.items():
        if '/' in t and key in [x.strip() for x in t.split('/')]:
            return v
    raise KeyError('no photoionization cross-section for level ' + levname)


def cross_section(term, nu, edge, nef, zion, hyd):
    """cross-section [cm^2] at frequency nu [1e15 Hz] for a level with
    ionization edge 'edge' [1e15 Hz] (SUB_PHOT_GEN)."""
    typ, a, unorm = term
    if typ == 0:
        return 0.0                        # no photoionization from this term
    if nu < edge / (1 + 1e-14):
        return 0.0
    u = nu / edge
    if typ in (20, 21):
        if u >= unorm[-1]:
            s = a[-1] * (unorm[-1] / u) ** 3
        else:
            j = np.searchsorted(unorm, u, side='right') - 1
            j = min(max(j, 0), len(unorm) - 2)
            t1 = (unorm[j + 1] - u) / (unorm[j + 1] - unorm[j])
            s = (1 - t1) * a[j + 1] + t1 * a[j]
        return max(s, 0.0) * MBARN
    if typ == 1:
        if a[0] == 0:
            return 0.0
        ru = edge / nu
        return a[0] * (a[1] + (1 - a[1]) * ru) * ru ** a[2] * MBARN
    if typ == 7:
        if a[0] == 0:
            return 0.0
        ru = (edge + a[3]) / nu
        if ru > 1:
            return 0.0
        return a[0] * (a[1] + (1 - a[1]) * ru) * ru ** a[2] * MBARN
    if typ == 4:
        ru = edge / nu
        s = a[0] + ru * (a[1] + ru * (a[2] + ru * (a[3] + ru * (a[4] + ru * a[5]))))
        return max(s, 0.0) * MBARN
    if typ == 5:
        x = math.log10(min(u, a[4]))
        s = 10 ** (a[0] + x * (a[1] + x * (a[2] + x * a[3])))
        if u > a[4]:
            s *= (a[4] / u) ** 2
        return s * MBARN
    if typ == 6:
        x = math.log10(u * (1 + 3e-16))
        if x < 0:
            return 0.0
        if x < a[4]:
            t1 = ((a[3] * x + a[2]) * x + a[1]) * x + a[0]
        else:
            t1 = a[5] + a[6] * x
        return 10 ** t1 * MBARN
    if typ == 9:
        s = 0.0
        e = edge
        for k in range(len(a) // 8):
            p = a[8 * k:8 * k + 8]
            if k > 0:
                e = EV_TO_HZ15 * p[2]
            if nu >= e:
                uu = nu / p[3] / EV_TO_HZ15
                t1 = (uu - 1) ** 2 + p[7] ** 2
                t2 = uu ** (5.5 + p[1] - 0.5 * p[6])
                t3 = (1 + math.sqrt(uu / p[5])) ** p[6]
                s += t1 * p[4] / t2 / t3
        return s * MBARN
    if typ == 3:
        n = int(round(a[1]))
        g = hyd.gaunt_n(n, u)
        alpha_bf = 2.815e-6 * zion ** 4
        s = alpha_bf * g * a[0] / nef / n / (nu * nef) ** 3   # program units
        return s * 1e-10
    if typ in (2, 8):
        n, lst, lend = int(round(a[0])), int(round(a[1])), int(round(a[2]))
        if typ == 8:
            e8 = edge + a[3]
            if nu < e8:
                return 0.0
            u = nu / e8
        tot = 0.0
        for l in range(lst, lend + 1):
            tot += (2 * l + 1) * 10 ** hyd.log_cross_l(n, l, u)
        if typ == 2:
            tot *= (nef / (n * zion)) ** 2
        else:
            tot /= zion ** 2
        tot /= (lend - lst + 1) * (lend + lst + 1)
        return tot * 1e-10
    raise NotImplementedError('cross-section type %d' % typ)


def u_grid(term):
    typ, a, unorm = term
    if typ == 0:
        return np.array([1.0, 1e3])
    if typ in (20, 21):
        u = np.array(unorm, dtype=float)
        u = u[u >= 1.0]
        if len(u) == 0 or u[0] > 1.0:
            u = np.concatenate([[1.0], u])
        return u
    # analytic: logarithmic grid to u = 1000, finer near the edge
    return np.unique(np.concatenate([np.geomspace(1.0, 2.0, 41), np.geomspace(2.0, 1e3, 120)]))


# ---------------------------------------------------------------------------
# collision file (GEN_OMEGA_RD_V2)
# ---------------------------------------------------------------------------
def read_col(path, levnames, gw):
    """Collision strengths (GEN_OMEGA_RD_V2).  A name without [J] matches
    the level of that name, or all its J components (multiplet collision
    strengths are split as Omega(i,j) = g_i g_j Omega / (G_L G_U)); 'I'
    denotes ionization (Omega(i,I) = g_i Omega / G_L)."""
    lines = read_lines(path)
    ntr = header_value(lines, '!Number of transitions', lambda s: int(float(s)))
    nt = header_value(lines, '!Number of T values OMEGA tabulated at', lambda s: int(float(s)))
    oscale = 1.0
    oset = 0.1
    for s in lines:
        if '!Scaling factor for OMEGA' in s:
            oscale = fnum(s.split('!')[0].split()[0])
        if '!Value for OMEGA if f=0' in s:
            oset = fnum(s.split('!')[0].split()[0])
    k = next(i for i, s in enumerate(lines) if 'Transition\\T' in s)
    tvals = [fnum(t) for t in lines[k].split('Transition\\T')[1].split()][:nt]
    assert len(tvals) == nt, (len(tvals), nt)
    exact = {nm: i + 1 for i, nm in enumerate(levnames)}
    base = {}
    for i, nm in enumerate(levnames):
        b = nm.split('[')[0]
        base.setdefault(b, []).append(i + 1)

    def match(name):
        kk = name.find('[')
        if kk >= 0:
            if name in exact:
                return [exact[name]]
            if name[:kk] in exact:
                return [exact[name[:kk]]]
            return []
        if name in exact:
            return [exact[name]]
        return list(base.get(name, []))

    om = {}
    omi = {}
    nread = 0
    for s in lines[k + 1:]:
        if not s.strip() or s.strip().startswith('!'):
            continue
        if '-' not in s:
            continue
        nread += 1
        if nread > ntr:
            break
        st = s.strip()
        j = st.index('-')
        low = st[:j].strip()
        rest = st[j + 1:].lstrip()
        jj = rest.find('  ')
        if jj < 0:
            continue
        up = rest[:jj].strip()
        try:
            vals = [fnum(x) for x in rest[jj:].split()[:nt]]
        except ValueError:
            continue
        if len(vals) < nt:
            continue
        lows = match(low)
        if not lows:
            continue
        gl = sum(gw[i - 1] for i in lows)
        if up == 'I':
            for i in lows:
                o = omi.setdefault(i, [0.0] * nt)
                for t in range(nt):
                    o[t] += vals[t] * gw[i - 1] / gl
            continue
        ups = match(up)
        if not ups:
            continue
        gu = sum(gw[i - 1] for i in ups)
        for i in lows:
            for u in ups:
                if i == u:
                    continue
                key = (min(i, u), max(i, u))
                o = om.setdefault(key, [0.0] * nt)
                for t in range(nt):
                    o[t] += vals[t] * gw[i - 1] * gw[u - 1] / (gl * gu)
    rows = [(a, b, v) for (a, b), v in sorted(om.items())]
    rows += [(a, 0, v) for a, v in sorted(omi.items())]
    return dict(t=tvals, scale=oscale, oset=oset, rows=rows, ntr=ntr)


# ---------------------------------------------------------------------------
# SAME_N (CMFGEN subs/same_n.f)
# ---------------------------------------------------------------------------
def posn_same(lev):
    """position (start, end; 1-based, 0 if not found) of the principal
    quantum number of the outer electron in a CMFGEN level name"""
    L = len(lev)

    def ch(i):                       # 1-based character, '' outside
        return lev[i - 1] if 1 <= i <= L else ''
    l1 = lev.find('__') + 1
    if l1 == 0:
        l1 = max(lev.find('SNG') + 1, lev.find('TRP') + 1)
    if l1 != 0:
        endn = l1 - 1
        stn = endn
        if stn > 1 and '0' < ch(stn - 1) < '9':
            stn -= 1
        if stn > 1 and '0' < ch(stn - 1) < '9':
            stn -= 1
        return stn, endn
    l1 = lev.find('[') + 1
    l1 = L if l1 == 0 else l1 - 1
    if ch(l1) in ('e', 'o'):
        l1 -= 1
    stn = endn = 0
    j = l1 - 3
    if ch(j + 1) < 'd' or ch(j + 1) > 's':
        j -= 1
    while stn == 0 and j > 0:
        if ch(j) == ')':
            while j > 0 and ch(j) != '(':
                j -= 1
            j -= 1
        if j > 0 and '0' <= ch(j) <= '9':
            if not (ch(j + 1) < 'd' or ch(j + 1) > 's'):
                endn = j
                stn = endn
                if stn > 1 and '0' <= ch(stn - 1) <= '9':
                    stn -= 1
                break
        j -= 1
    return stn, endn


def same_n(a, b):
    s1, e1 = posn_same(a)
    if s1 == 0:
        return False
    s2, e2 = posn_same(b)
    if s2 == 0:
        return False
    return a[s1 - 1:e1] == b[s2 - 1:e2] and a[:s1] == b[:s2]


# ---------------------------------------------------------------------------
TWO_PHOTON = {   # (element, ion stage): (lower name, upper name, A [1/s], g of decaying state)
    (1, 1): (['1___', '2s_2Se'], ['1___'], 8.225, 2.0),
    (2, 1): (['1s_2s_1Se', '1s2s1Se'], ['1s2_1Se', '1s1s1Se'], 51.9, 1.0),
    (2, 2): (['2___', '2s2Se'], ['1___', '1s2Se'], 526.4, 2.0),
}

ELEM = {1: 'H', 2: 'He', 6: 'C', 7: 'N', 8: 'O', 10: 'Ne', 11: 'Na', 12: 'Mg', 13: 'Al',
        14: 'Si', 16: 'S', 18: 'Ar', 20: 'Ca', 26: 'Fe', 27: 'Co', 28: 'Ni'}
ROMAN = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII']


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--z', type=int, required=True, help='atomic number')
    ap.add_argument('--ion', type=int, required=True, help='ion stage (1 = neutral)')
    ap.add_argument('--osc', required=True)
    ap.add_argument('--phot', required=True)
    ap.add_argument('--phot-b', action='append', default=[])
    ap.add_argument('--next-osc', default=None)
    ap.add_argument('--col', required=True)
    ap.add_argument('--hyd', default=None)
    ap.add_argument('--ftos', default=None)
    ap.add_argument('--nlev', type=int, default=0, help='levels to keep (0: all)')
    ap.add_argument('-o', '--out', required=True)
    args = ap.parse_args()

    osc = read_osc(args.osc)
    nlev = osc['nlev'] if args.nlev <= 0 else min(args.nlev, osc['nlev'])
    levs = osc['levels'][:nlev]
    names = [lv['name'] for lv in levs]
    chi = osc['chi']
    zion = osc['zion']
    sl = read_ftos(args.ftos, nlev) if args.ftos else list(range(1, nlev + 1))

    hdr, terms = read_phot(args.phot)
    split_j = hdr.get('!Split J levels', 'False').strip().lower().startswith('t')
    gion = fnum(hdr.get('!Statistical weight of ion', '1.0'))
    exc = fnum(hdr.get('!Excitation energy of final state', '0.0'))
    # route A final state above the ground state of the next ion [cm^-1]:
    # the cross-sections are functions of nu/(edge + exc)
    exc15 = exc * CLIGHT * 1e-15
    nroutes = int(fnum(hdr.get('!Number of photoionization routes', '1')))
    # further routes: (terms, split_j, excitation of the final term [1e15 Hz])
    routes = []
    nxt = read_osc(args.next_osc)['levels'] if args.next_osc else None
    for pb in args.phot_b:
        hb, tb = read_phot(pb)
        fin = None
        for s0 in read_lines(pb):
            if '!Final state in ion' in s0:
                fin = s0.split('!')[0].strip()
                break
        eb = final_term_energy(fin, nxt) if (nxt is not None and fin) else None
        if eb is None:
            print('note: final state %s of %s not found in the next ion: route skipped' % (fin, pb),
                  file=sys.stderr)
            continue
        routes.append((tb, hb.get('!Split J levels', 'False').strip().lower().startswith('t'),
                       eb * CLIGHT * 1e-15, fin))
        print('note: route %s to %s at %.2f cm^-1 added' % (os.path.basename(pb), fin, eb), file=sys.stderr)
    if nroutes > 1 + len(routes):
        print('note: %d photoionization routes in the data, %d used' % (nroutes, 1 + len(routes)),
              file=sys.stderr)
    hyd = Hyd(args.hyd)
    nu_inf = 1e-15 * RINF * CLIGHT
    nu_inf /= (1 + 5.48597e-4) if args.z == 1 else (1 + 5.48597e-4 / (2 * args.z))

    trans = [t for t in osc['trans'] if t['u'] <= nlev and t['l'] < t['u']]
    col = read_col(args.col, names, [lv['g'] for lv in levs])

    ename = ELEM.get(args.z, 'Z%d' % args.z)
    with open(args.out, 'w') as f:
        w = f.write
        w('# SuperLite NLTE model atom, format 1\n')
        w('# %s %s converted from CMFGEN atomic data on %s\n' % (
            ename, ROMAN[args.ion - 1], datetime.date.today().isoformat()))
        w('#   osc:  %s\n#   phot: %s\n#   col:  %s\n' % tuple(
            os.path.relpath(p, os.path.dirname(os.path.dirname(os.path.dirname(p))))
            for p in (args.osc, args.phot, args.col)))
        if args.ftos:
            w('#   f_to_s: %s\n' % os.path.basename(args.ftos))
        w('element %d %s\n' % (args.z, ename))
        w('ion %d\n' % args.ion)
        w('zion %.4f\n' % zion)
        w('chi %.6f            ! ionization energy of level 1 [cm^-1]\n' % chi)
        w('gion %.4f           ! statistical weight of the ground state of the next ion\n' % gion)
        w('nlev %d\n' % nlev)
        w('# level name g E[cm^-1] superlevel\n')
        for i, lv in enumerate(levs):
            w('%5d %-24s %8.2f %16.6f %5d\n' % (i + 1, lv['name'], lv['g'], lv['e'], sl[i]))
        w('nline %d\n' % len(trans))
        w('# lower upper f A[1/s] lambda[A](vac<2000<air; <0: not measured) same_n\n')
        for t in trans:
            sn = 1 if same_n(names[t['l'] - 1], names[t['u'] - 1]) else 0
            w('%5d %5d %12.4e %12.4e %14.4f %d\n' % (t['l'], t['u'], t['f'], t['a'], t['lam'], sn))
        # two-photon
        tp = TWO_PHOTON.get((args.z, args.ion))
        tpl = []
        if tp:
            up = next((i + 1 for i, nm in enumerate(names) if nm in tp[0]), 0)
            lo = next((i + 1 for i, nm in enumerate(names) if nm in tp[1]), 0)
            if up and lo:
                # rate per atom of the (possibly l-bundled) level
                tpl.append((lo, up, tp[2] * tp[3] / levs[up - 1]['g']))
        w('ntwo %d\n' % len(tpl))
        for lo, up, a in tpl:
            w('%5d %5d %12.4e\n' % (lo, up, a))
        # collisions
        rows = [r for r in col['rows'] if r[0] <= nlev and r[1] <= nlev]
        w('ncolt %d\n' % len(col['t']))
        w('colt ' + ' '.join('%.4g' % t for t in col['t']) + '   ! [10^4 K]\n')
        w('omegaset %.4g\n' % col['oset'])
        w('omegascale %.4g\n' % col['scale'])
        w('ncol %d\n' % len(rows))
        w('# lower upper(0: ionization) Omega(T)\n')
        for lo, up, om in rows:
            w('%5d %5d ' % (lo, up) + ' '.join('%11.4e' % x for x in om) + '\n')
        # photoionization
        w('nphot %d\n' % nlev)
        w('# level npts, then u=nu/nu_edge sigma[cm^2]\n')
        for i, lv in enumerate(levs):
            edge0 = (chi - lv['e']) * CLIGHT * 1e-15
            edge = edge0 + exc15
            nef = zion * math.sqrt(nu_inf / edge) if edge > 0 else 0.0
            term = term_for_level(lv['name'], terms, split_j)
            ug = u_grid(term)
            sig = np.array([cross_section(term, u * edge, edge, nef, zion, hyd) for u in ug])
            ug = ug * edge / edge0              # tabulated against nu/nu0 of the level
            if sig[0] <= 0 and len(sig) > 1:
                # threshold value as the first tabulated value above it (CMFGEN SET_TO_EDGE)
                sig[0] = sig[1]
            # further routes: summed on the union of the grids
            extra = []
            for tb, sjb, eb15, fin in routes:
                try:
                    termb = term_for_level(lv['name'], tb, sjb)
                except KeyError:
                    continue
                if termb[0] == 0:
                    continue                    # no cross-section on this route
                edgeb = edge0 + eb15
                nefb = zion * math.sqrt(nu_inf / edgeb)
                ugb = u_grid(termb)
                sgb = np.array([cross_section(termb, u * edgeb, edgeb, nefb, zion, hyd) for u in ugb])
                if sgb[0] <= 0 and len(sgb) > 1:
                    sgb[0] = sgb[1]
                extra.append((termb, edgeb, nefb, ugb * edgeb / edge0, sgb))
            if extra:
                grid = [ug] + [e[3] for e in extra] + [np.array([e[3][0] * (1 - 1e-7) for e in extra])]
                # the table is interpolated linearly in u: where a route is in
                # its nu^-3 tail, 40 points per decade up to the last point
                ends = [ug[-1]] + [e[3][-1] for e in extra]
                if max(ends) > 1.05 * min(ends):
                    grid.append(np.geomspace(min(ends), max(ends),
                                             int(40 * math.log10(max(ends) / min(ends))) + 2))
                ugm = np.unique(np.concatenate(grid))
                sigm = np.interp(ugm, ug, sig, right=np.nan)
                tail = ugm > ug[-1]
                sigm[tail] = sig[-1] * (ug[-1] / ugm[tail]) ** 3
                for termb, edgeb, nefb, ub, sb in extra:
                    add = np.interp(ugm, ub, sb, left=0.0, right=np.nan)
                    tl = ugm > ub[-1]
                    add[tl] = sb[-1] * (ub[-1] / ugm[tl]) ** 3
                    add[ugm < ub[0]] = 0.0
                    sigm = sigm + add
                ug, sig = ugm, sigm
            w('%5d %5d\n' % (i + 1, len(ug)))
            for u, s in zip(ug, sig):
                w('%14.7e %12.5e\n' % (u, s))
    print('wrote %s: %d levels, %d lines, %d collision rows' % (args.out, nlev, len(trans), len(rows)))


if __name__ == '__main__':
    main()
