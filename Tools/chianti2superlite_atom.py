#!/usr/bin/env python3
# This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
"""SuperLite NLTE model atom (format 1, as Tools/cmfgen2superlite_atom.py
writes) of an ion from the CHIANTI database: for the highly charged ions
that CMFGEN's atomic data do not have (Fe XVII-XXVI, the H- and He-like
ions of Ne-Ca, ...).

  levels       .elvlc: observed energies (theoretical where none), g = 2J+1;
               the levels below the ionization limit, lowest first (at most
               --nlev); super-levels: the levels of one term (configuration,
               2S+1, L), or none (--no-super)
  lines        .wgfa: f from A (A > 0; f_lu = m_e c / (8 pi^2 e^2) lambda^2
               g_u/g_l A), the frequency from the level energies; those with
               g_l f >= --gfmin or from the lowest --gflev levels (CMFGEN's
               GF_CUT and GF_LEV_CUT); the
               zero-wavelength decays of the ground state's neighbours
               (H-like 2s, He-like 1s2s 1S0) as two-photon decays
  collisions   .scups: effective collision strengths (Burgess & Tully 1992
               scaling, as CHIANTI's descale) on a temperature grid; pairs
               without them get SuperLite's rules (van Regemorter, omegaset)
  photoioniz.  ground state: the valence shells of Verner & Yakovlev (1995)
               (Data/data.bf_verner; the core shells, 1s for >= 3 electrons,
               1s-2p for >= 11, 1s-3p for >= 19, are SuperLite's Auger
               absorption, in_nltex_auger); excited levels: Kramers,
               sigma = 7.91e-18 n*/z^2 (nu_0/nu)^3 cm^2, n* the effective
               quantum number of the level's edge
  dielectronic the total dielectronic recombination of the next ion
  recombin.    (its .drparams, Badnell's fit T^-1.5 sum_i c_i exp(-E_i/T)) as
               resonances in the ground-state cross-section: by the Milne
               relation a resonance at h nu_0 + k E_i with the integrated
               cross-section c_i / [g_1/(2 g_+) C_S 8 pi nu_i^2/c^2] (C_S =
               (h^2/(2 pi m_e k))^1.5) gives exactly the term i of the fit
               (as CMFGEN's dielectronic resonances); Gaussians of --vsm km/s,
               narrower close to the edge (all of a resonance above it)

The radiative recombination that the Milne relation gives from these
cross-sections is scaled to CHIANTI's total (the next ion's .rrparams,
Badnell 2006, all n) at 3e4 K through the excited levels' cross-sections
(the factor in the atom's header) and printed against it.

usage: chianti2superlite_atom.py <CHIANTI dbase dir> <element symbol> <stage>
           <out.atom> [--nlev 150] [--no-super] [--vsm 1000] [--gfmin 1e-3] [--gflev 9]
e.g.   chianti2superlite_atom.py $XUVTOP fe 17 fe17.atom
The atoms are made from the CHIANTI database (Dere et al. 1997; Del Zanna
et al. 2021): cite CHIANTI when using them; they are not distributed with
SuperLite.
"""
import argparse
import datetime
import math
import os
import re
import sys

import numpy as np
from scipy.interpolate import splev, splrep

CLIGHT = 2.99792458e10
HH = 6.6260693e-27
KB = 1.3806505e-16
ME = 9.1093826e-28
EE = 4.80320441e-10
RYD_CM = 109737.31568
RYD_K = 157887.51
SAHAC = 4.1412e-16                    # (h^2/(2 pi m_e k))^1.5 [cm^3 K^1.5]
ROMAN = ['I', 'II', 'III', 'IV', 'V', 'VI', 'VII', 'VIII', 'IX', 'X', 'XI', 'XII', 'XIII',
         'XIV', 'XV', 'XVI', 'XVII', 'XVIII', 'XIX', 'XX', 'XXI', 'XXII', 'XXIII', 'XXIV',
         'XXV', 'XXVI', 'XXVII', 'XXVIII', 'XXIX', 'XXX']
SYM = ['h', 'he', 'li', 'be', 'b', 'c', 'n', 'o', 'f', 'ne', 'na', 'mg', 'al', 'si', 'p',
       's', 'cl', 'ar', 'k', 'ca', 'sc', 'ti', 'v', 'cr', 'mn', 'fe', 'co', 'ni', 'cu', 'zn']
LETTER_L = 'SPDFGHIKLMNOQRTUVWXYZ'
COLT = [0.3, 0.5, 1, 2, 3, 5, 10, 20, 30, 50, 100, 200, 300, 500, 1000]   # [1e4 K]
VERNER_L = {1: 0, 2: 0, 3: 1, 4: 0, 5: 1, 6: 2, 7: 0}                      # shell -> l


def records(path):
    """data lines of a CHIANTI file up to its '-1' terminator"""
    out = []
    for s in open(path):
        if s.strip().startswith('-1'):
            break
        out.append(s.rstrip('\n'))
    return out


def read_elvlc(path):
    levs = []
    num = r'(-?\d+\.?\d*(?:[eE][-+]?\d+)?)'
    pat = re.compile(r'^\s*(\d+)\s+(.*?)\s+(\d+)\s+([A-Z])\s+(\d+\.?\d*)\s+' + num + r'\s+' + num)
    for s in records(path):
        m = pat.match(s)                     # (files may add columns at the end)
        if not m:
            continue
        idx = int(m.group(1))
        conf = m.group(2).split()
        if len(conf) > 1 and conf[-1].isdigit():
            conf = conf[:-1]                 # (a numeric label before 2S+1)
        conf = '.'.join(conf)
        mult, l2s, j = int(m.group(3)), m.group(4), float(m.group(5))
        eobs, eth = float(m.group(6)), float(m.group(7))
        e = eobs if (eobs > 0 or idx == 1) else eth
        levs.append(dict(idx=idx, conf=conf, mult=mult, L=l2s, j=j, g=2 * j + 1, e=e,
                         name='%s_%d%s[%g]' % (conf, mult, l2s, j)))
    return levs


def read_wgfa(path):
    out = []
    for s in records(path):
        w = s.split()
        if len(w) < 5:
            continue
        try:
            out.append((int(w[0]), int(w[1]), float(w[2]), float(w[3]), float(w[4])))
        except ValueError:
            continue
    return out


def read_scups(path):
    """(l, u, de [Ryd], ttype, c, x, y) per transition"""
    lines = records(path)
    out = []
    i = 0
    while i + 2 < len(lines):
        w = lines[i].split()
        if len(w) < 8:
            i += 1
            continue
        l, u = int(w[0]), int(w[1])
        de, nspl, ttype, c = float(w[2]), int(w[5]), int(w[6]), float(w[7])
        x = np.array([float(v) for v in lines[i + 1].split()][:nspl])
        y = np.array([float(v) for v in lines[i + 2].split()][:nspl])
        out.append((l, u, de, ttype, c, x, y))
        i += 3
    return out


def upsilon(de, ttype, c, x, y, temp):
    """Burgess & Tully descaling (CHIANTI descale_bt) at temperatures temp [K]"""
    kte = temp / RYD_K / de
    if ttype in (1, 4):
        xt = 1 - np.log(c) / np.log(kte + c)
    else:
        xt = kte / (kte + c)
    tck = splrep(x, y, s=0, k=min(3, len(x) - 1))
    yt = splev(np.clip(xt, x[0], x[-1]), tck)
    if ttype == 1:
        return yt * np.log(kte + math.e)
    if ttype == 2:
        return yt
    if ttype == 3:
        return yt / (kte + 1)
    if ttype == 4:
        return yt * np.log(kte + c)
    if ttype == 5:
        return yt / kte
    if ttype == 6:
        return 10 ** yt
    raise ValueError('scups type %d' % ttype)


def read_ip(dbase, z, stage):
    for s in open(os.path.join(dbase, 'ip', 'chianti.ip')):
        w = s.split()
        if len(w) == 3 and w[0].lstrip('-').isdigit() and int(w[0]) == z and int(w[1]) == stage:
            return float(w[2])
    raise KeyError('no ionization energy for Z=%d stage %d' % (z, stage))


def read_drparams(path):
    if not os.path.isfile(path):
        return None
    lines = records(path)
    if int(lines[0].split()[0]) != 1:
        raise ValueError('%s: only Badnell fits (type 1)' % path)
    e = [float(v) for v in lines[1].split()[2:]]
    c = [float(v) for v in lines[2].split()[2:]]
    return [(ei, ci) for ei, ci in zip(e, c) if ci > 0]


def rr_chianti(path, temp):
    lines = records(path)
    typ = int(lines[0].split()[0])
    w = [float(v) for v in lines[1].split()]
    if typ in (1, 2):
        a, b, t0, t1 = w[3:7]
        if typ == 2:
            b = b + w[7] * np.exp(-w[8] / temp)
        return a / (np.sqrt(temp / t0) * (1 + np.sqrt(temp / t0)) ** (1 - b)
                    * (1 + np.sqrt(temp / t1)) ** (1 + b))
    if typ == 3:
        return w[2] * (temp / 1e4) ** (-w[3])
    return np.full_like(temp, np.nan)


def verner_valence(path, z, nelec):
    """the valence shells (Verner & Yakovlev 1995) of the ion: E_thr, E_0,
    sigma_0, y_a, P, y_w, l"""
    ninn = 0 if nelec < 3 else 1 if nelec < 11 else 3 if nelec < 19 else 5
    out = []
    for s in open(path):
        if s.startswith('# Verner 1996'):
            break
        if s.startswith('#'):
            continue
        w = s.split()
        ns, ne_, zz = int(w[0]), int(w[1]), int(w[2])
        if zz == z and ne_ == nelec and ns > ninn:
            out.append([float(v) for v in w[3:9]] + [VERNER_L[ns]])
    return out


def sigma_verner(shells, e_ev):
    s = np.zeros_like(e_ev)
    for ethr, e0, s0, ya, p, yw, l in shells:
        y = e_ev / e0
        f = ((y - 1) ** 2 + yw ** 2) * y ** (0.5 * p - 5.5 - l) * (1 + np.sqrt(y / ya)) ** (-p)
        s += np.where(e_ev >= ethr * (1 - 1e-6), s0 * f, 0.0)
    return s * 1e-18


TNORM = 3e4


def milne(phot, nu0, g, gion, shells, ip, dr, t, full=False):
    """recombination [cm^3/s] at t by the Milne relation: radiative into
    the ground state, into the excited levels, and (full) the dielectronic
    part (the resonances of the ground state)"""
    ev_cm = 8065.544
    out = [0.0, 0.0, 0.0]
    for k, (ub, s) in enumerate(phot):
        s0 = s
        if k == 0 and dr:
            s0 = sigma_verner(shells, np.maximum(ub * ip / ev_cm, shells[0][0]))
        xk = np.linspace(0, 60, 3001) * KB * t / (HH * nu0[k])
        uf = np.unique(np.concatenate([1 + xk[1 + xk <= 1000], ub]))
        nu = uf * nu0[k]
        w = 8 * math.pi * nu ** 2 / CLIGHT ** 2 * np.exp(-HH * (nu - nu0[k]) / (KB * t))
        fac = g[k] / (2 * gion) * SAHAC * t ** -1.5
        sf0 = np.exp(np.interp(np.log(uf), np.log(ub), np.log(np.maximum(s0, 1e-99))))
        r = fac * np.trapezoid(sf0 * w, nu)
        out[0 if k == 0 else 1] += r
        if full and k == 0 and dr:
            out[2] += fac * np.trapezoid((np.interp(uf, ub, s) - np.interp(uf, ub, s0)) * w, nu)
    return tuple(out) if full else (out[0], out[1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('dbase')
    ap.add_argument('element')
    ap.add_argument('stage', type=int)
    ap.add_argument('out')
    ap.add_argument('--nlev', type=int, default=150)
    ap.add_argument('--no-super', action='store_true')
    ap.add_argument('--vsm', type=float, default=1000.)
    ap.add_argument('--gfmin', type=float, default=1e-3)
    ap.add_argument('--gflev', type=int, default=9)
    ap.add_argument('--verner', default=os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                                     '..', 'Data', 'data.bf_verner'))
    a = ap.parse_args()
    el = a.element.lower()
    z = SYM.index(el) + 1
    st = a.stage
    nelec = z - st + 1
    name = '%s_%d' % (el, st)
    d = os.path.join(a.dbase, el, name)
    ip = read_ip(a.dbase, z, st)                                   # [cm^-1]
    # levels below the limit, lowest first
    allv = read_elvlc(os.path.join(d, name + '.elvlc'))
    allv = sorted([v for v in allv if v['e'] < ip], key=lambda v: (v['e'], v['idx']))[:a.nlev]
    keep = {v['idx']: i for i, v in enumerate(allv)}
    n = len(allv)
    # super-levels: the terms in order of their lowest level
    sl, terms = [], {}
    for v in allv:
        key = (v['conf'], v['mult'], v['L']) if not a.no_super else v['idx']
        if key not in terms:
            terms[key] = len(terms) + 1
        sl.append(terms[key])
    # next ion
    if st == z:
        gion = 1.0
    else:
        nxt = read_elvlc(os.path.join(a.dbase, el, '%s_%d' % (el, st + 1), '%s_%d.elvlc' % (el, st + 1)))
        n1 = min(nxt, key=lambda v: v['e'])
        # (the ground term, as CMFGEN's cross-sections lead to it)
        gion = sum(v['g'] for v in nxt if (v['conf'], v['mult'], v['L']) == (n1['conf'], n1['mult'], n1['L']))
    g = np.array([v['g'] for v in allv])
    ecm = np.array([v['e'] for v in allv]) - allv[0]['e']
    # lines and two-photon decays
    lines, two = [], []
    for i, j, wvl, gf, aval in read_wgfa(os.path.join(d, name + '.wgfa')):
        if i not in keep or j not in keep or aval <= 0:
            continue
        l, u = keep[i], keep[j]
        if ecm[u] < ecm[l]:
            l, u = u, l
        if wvl == 0:
            if l == 0 and u <= 2:
                two.append((l + 1, u + 1, aval))
            continue
        if ecm[u] <= ecm[l]:
            continue
        lam = 1.0 / (ecm[u] - ecm[l])                                  # [cm]
        f = ME * CLIGHT / (8 * math.pi ** 2 * EE ** 2) * lam ** 2 * g[u] / g[l] * aval
        if g[l] * f < a.gfmin and l + 1 > a.gflev:
            continue
        lines.append((l + 1, u + 1, f, aval, lam * 1e8))
    lines.sort()
    # collisions
    tgrid = np.array(COLT) * 1e4
    cols = []
    scf = os.path.join(d, name + '.scups')
    if os.path.isfile(scf):
        for l, u, de, ttype, c, x, y in read_scups(scf):
            if l not in keep or u not in keep or de <= 0:
                continue
            lo, up = sorted((keep[l], keep[u]))
            if lo == up:
                continue
            ups = upsilon(de, ttype, c, x, y, tgrid)
            cols.append((lo + 1, up + 1, np.maximum(ups, 0.0)))
    cols.sort(key=lambda r: (r[0], r[1]))
    # photoionization: grid in u = nu/nu_edge to 1000
    shells = verner_valence(a.verner, z, nelec)
    if not shells:
        sys.exit('no Verner valence shells for Z=%d, %d electrons' % (z, nelec))
    ev_cm = 8065.544
    dr = read_drparams(os.path.join(a.dbase, el, '%s_%d' % (el, st + 1), '%s_%d.drparams' % (el, st + 1))) \
        if st < z else None
    nu0 = (ip - ecm) * CLIGHT                                          # [Hz]
    phot = []
    for k in range(n):
        ub = np.unique(np.concatenate([np.geomspace(1, 1000, 120), [1.0 + 1e-6]]))
        if k == 0:
            res = []
            if dr:
                for ek, ck in dr:
                    nur = nu0[0] + KB * ek / HH
                    s_int = ck / (g[0] / (2 * gion) * SAHAC * 8 * math.pi * nur ** 2 / CLIGHT ** 2)
                    res.append((nur / nu0[0], s_int))
            # (width: --vsm, or a sixth of the distance from the edge, so that
            # a resonance close to the edge stays above it)
            dl0 = a.vsm / 2.99792458e5
            res = [(u0, si, max(min(dl0, (u0 - 1) / (6 * u0)), 1e-6)) for u0, si in res]
            nodes = [ub] + [u0 * (1 + dl * np.linspace(-5, 5, 41)) for u0, _, dl in res if u0 < 1000]
            ub = np.unique(np.concatenate(nodes))
            ub = ub[(ub >= 1) & (ub <= 1000)]
            s = sigma_verner(shells, ub * ip / ev_cm)
            if s[0] <= 0:                                   # (threshold above the CHIANTI limit)
                s = sigma_verner(shells, np.maximum(ub * ip / ev_cm, shells[0][0]))
            for u0, s_int, dl in res:
                if u0 >= 1000:
                    continue
                x = (ub - u0) / (u0 * dl)
                phi = np.where(np.abs(x) <= 5, np.exp(-x * x), 0.0)
                area = np.trapezoid(phi, ub) * nu0[0]
                if area > 0:
                    s = s + s_int / area * phi
        else:
            nstar = st * math.sqrt(RYD_CM / (ip - ecm[k]))
            s = 7.91e-18 * nstar / st ** 2 * ub ** -3.0
        phot.append((ub, s))
    # Milne radiative (and dielectronic) recombination against CHIANTI; the
    # excited levels' cross-sections scaled so that the radiative
    # recombination matches CHIANTI's total at TNORM (within 0.05-20)
    nname = '%s_%d' % (el, st + 1)
    rrf = os.path.join(a.dbase, el, nname, nname + '.rrparams')       # (X+1 recombining to X)
    scale = 1.0
    if os.path.isfile(rrf) and n > 1:
        a0, ae = milne(phot, nu0, g, gion, shells, ip, dr, TNORM)
        target = rr_chianti(rrf, np.array([TNORM]))[0]
        if ae > 0 and target > a0:
            scale = min(max((target - a0) / ae, 0.05), 20.0)
        for k in range(1, n):
            phot[k] = (phot[k][0], phot[k][1] * scale)
    tt = np.array([1e4, 1e5, 1e6, 1e7])
    alpha = np.zeros_like(tt)
    alpha_rr = np.zeros_like(tt)
    for it, t in enumerate(tt):
        r0, re_, rd = milne(phot, nu0, g, gion, shells, ip, dr, t, full=True)
        alpha_rr[it] = r0 + re_
        alpha[it] = r0 + re_ + rd
    rrc = rr_chianti(rrf, tt) if os.path.isfile(rrf) else np.full_like(tt, np.nan)
    drc = np.array([t ** -1.5 * sum(c * math.exp(-e / t) for e, c in dr) for t in tt]) if dr else np.zeros_like(tt)
    # write
    with open(a.out, 'w') as f:
        w = f.write
        w('# SuperLite NLTE model atom, format 1\n')
        w('# %s %s from the CHIANTI database (v%s) on %s (Tools/chianti2superlite_atom.py)\n'
          % (el.capitalize(), ROMAN[st - 1], open(os.path.join(a.dbase, 'VERSION')).read().strip(),
             datetime.date.today().isoformat()))
        w('#   levels: %d of %d below the limit; super-levels: %s\n'
          % (n, len(read_elvlc(os.path.join(d, name + '.elvlc'))), 'none' if a.no_super else 'terms'))
        w('#   photoionization: Verner & Yakovlev 1995 valence shells (ground), Kramers (excited);'
          ' dielectronic: %s\n' % ('Badnell fit of %s_%d.drparams as resonances' % (el, st + 1) if dr else 'none'))
        w('#   excited levels: Kramers x %.3g (radiative recombination = CHIANTI .rrparams at %g K)\n'
          % (scale, TNORM))
        w('#   recombination [cm^3/s] at T = 1e4 1e5 1e6 1e7 K: radiative (Milne) %s; CHIANTI %s;'
          ' dielectronic (Milne) %s; CHIANTI %s\n'
          % (' '.join('%.2e' % v for v in alpha_rr), ' '.join('%.2e' % v for v in rrc),
             ' '.join('%.2e' % v for v in alpha - alpha_rr), ' '.join('%.2e' % v for v in drc)))
        w('element %d %s\n' % (z, el.capitalize()))
        w('ion %d\n' % st)
        w('zion %.4f\n' % st)
        w('chi %.6f            ! ionization energy of level 1 [cm^-1]\n' % (ip - allv[0]['e']))
        w('gion %.4f           ! statistical weight of the ground state of the next ion\n' % gion)
        w('nlev %d\n' % n)
        w('# level name g E[cm^-1] superlevel\n')
        for i, v in enumerate(allv):
            w('%5d %-24s %8.2f %16.6f %5d\n' % (i + 1, v['name'], v['g'], ecm[i], sl[i]))
        w('nline %d\n' % len(lines))
        w('# lower upper f A[1/s] lambda[A](vac<2000<air; <0: not measured) same_n\n')
        for l, u, fval, aval, lam in lines:
            w('%5d %5d %12.4e %12.4e %14.4f 0\n' % (l, u, fval, aval, lam))
        w('ntwo %d\n' % len(two))
        for l, u, aval in two:
            w('%5d %5d %12.4e\n' % (l, u, aval))
        w('ncolt %d\n' % len(COLT))
        w('colt ' + ' '.join('%g' % t for t in COLT) + '   ! [10^4 K]\n')
        w('omegaset 0.1\n')
        w('omegascale 1\n')
        w('ncol %d\n' % len(cols))
        w('# lower upper(0: ionization) Omega(T)\n')
        for l, u, ups in cols:
            w('%5d %5d ' % (l, u) + ' '.join('%11.4e' % v for v in ups) + '\n')
        w('nphot %d\n' % n)
        w('# level npts, then u=nu/nu_edge sigma[cm^2]\n')
        for k, (ub, s) in enumerate(phot):
            w('%5d %5d\n' % (k + 1, len(ub)))
            for uu, ss in zip(ub, s):
                w('%15.7e %12.5e\n' % (uu, ss))
    print('wrote %s: %d levels in %d super-levels, %d lines, %d two-photon, %d collision rows; '
          'recombination at 1e4/1e5/1e6 K: radiative %s (CHIANTI %s), dielectronic %s (CHIANTI %s)'
          % (a.out, n, max(sl), len(lines), len(two), len(cols),
             '/'.join('%.1e' % v for v in alpha_rr[:3]), '/'.join('%.1e' % v for v in rrc[:3]),
             '/'.join('%.1e' % v for v in (alpha - alpha_rr)[:3]), '/'.join('%.1e' % v for v in drc[:3])))


if __name__ == '__main__':
    main()
