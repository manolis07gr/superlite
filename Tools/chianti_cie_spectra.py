#!/usr/bin/env python3
# This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
"""
Emission spectra of hot gas in collisional ionization equilibrium (CIE) for
the radiative shocks of SuperLite (in_shock_cie): per element, the line,
free-free, free-bound and two-photon emissivity per n_e n_el
[erg cm^3 s^-1] in logarithmic wavelength bins, on a grid of temperatures,
from the CHIANTI database (Dere et al. 1997; Del Zanna et al. 2021) with
ChiantiPy (CIE ion fractions of CHIANTI's default ioneq file).

usage (python with ChiantiPy, XUVTOP pointing to the CHIANTI database):
  python3 -I chianti_cie_spectra.py -o shock_cie.dat [--elements h,he,c,...]
          [--logt 5.0,8.0,0.1] [--ne 1e10] [--wl 1,3e5,0.004] [--fmin 1e-6]

The table is not distributed with SuperLite (see the CHIANTI terms of use);
cite CHIANTI when using it.  An ion contributes at the temperatures where
its CIE fraction exceeds fmin.

Output (text): header lines starting with '#', then
  nel <n>        elements <Z ...>
  nt <n>         logt <values>
  nw <n>         edges [A] (nw+1 values)
  ne <n_e>
and, per element and temperature, a line 'Z logT' followed by the nw bin
emissivities 4 pi eps [erg cm^3 s^-1] per n_e n_el.
"""
import argparse
import os
import sys
import time
import warnings

import numpy as np

ZOF = {'h': 1, 'he': 2, 'c': 6, 'n': 7, 'o': 8, 'ne': 10, 'na': 11, 'mg': 12, 'al': 13, 'si': 14,
       's': 16, 'ar': 18, 'ca': 20, 'fe': 26, 'ni': 28}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--elements', default='h,he,c,n,o,ne,mg,si,s,ar,ca,fe,ni')
    ap.add_argument('--logt', default='5.0,8.0,0.1')
    ap.add_argument('--ne', type=float, default=1e10)
    ap.add_argument('--wl', default='1,3e5,0.004', help='min, max [A], bin width [dex]')
    ap.add_argument('--fmin', type=float, default=1e-6)
    a = ap.parse_args()
    if 'XUVTOP' not in os.environ:
        sys.exit('set XUVTOP to the CHIANTI database directory')
    warnings.simplefilter('ignore')
    import ChiantiPy.core as ch
    import ChiantiPy.tools.io as chio
    import ChiantiPy.tools.data as chdata

    t0, t1, dt = (float(x) for x in a.logt.split(','))
    logt = np.round(np.arange(t0, t1 + 0.5 * dt, dt), 4)
    temp = 10.0 ** logt
    w0, w1, dw = (float(x) for x in a.wl.split(','))
    edges = 10.0 ** np.arange(np.log10(w0), np.log10(w1) + 0.5 * dw, dw)
    nw = len(edges) - 1
    wmid = np.sqrt(edges[:-1] * edges[1:])
    dwl = np.diff(edges)
    els = a.elements.split(',')
    ioneq_all = chio.ioneqRead(ioneqName=chdata.Defaults['ioneqfile'])   # default CIE fractions
    print('ion fractions:', chdata.Defaults['ioneqfile'], file=sys.stderr)
    tab = {}
    for el in els:
        z = ZOF[el]
        em = np.zeros((len(temp), nw))
        tstart = time.time()
        for stage in range(1, z + 2):
            ions = '%s_%d' % (el, stage)
            frac = np.interp(logt, np.log10(ioneq_all['ioneqTemperature']),
                             ioneq_all['ioneqAll'][z - 1, stage - 1], left=0.0, right=0.0)
            sel = frac > a.fmin
            if not sel.any():
                continue
            ts = temp[sel]
            # lines of the ion (not the bare nucleus)
            if stage <= z:
                try:
                    if ts.size == 1:     # (ChiantiPy fails for one-element arrays)
                        ion = ch.ion(ions, temperature=float(ts[0]), eDensity=a.ne, abundance=1.0)
                    else:
                        ion = ch.ion(ions, temperature=ts, eDensity=np.full(ts.size, a.ne), abundance=1.0)
                    ion.intensity()
                    it = ion.Intensity
                    if 'errorMessage' not in it:
                        inten = np.atleast_2d(it['intensity'])          # (ntemp, nlines)
                        wv = np.abs(np.asarray(it['wvl'], float))
                        k = np.searchsorted(edges, wv) - 1
                        ok = (k >= 0) & (k < nw)
                        for j, ij in enumerate(np.where(sel)[0]):
                            em[ij] += 4 * np.pi * np.bincount(k[ok], weights=inten[j][ok], minlength=nw)
                    # two-photon continuum (H- and He-like ions)
                    if hasattr(ion, 'twoPhoton'):
                        try:
                            ion.twoPhoton(wmid)
                            tp = np.atleast_2d(ion.TwoPhoton['intensity'])
                            for j, ij in enumerate(np.where(sel)[0]):
                                em[ij] += 4 * np.pi * tp[j] * dwl
                        except Exception:
                            pass
                except Exception as e:
                    print('  %s lines: %s' % (ions, str(e)[:80]), file=sys.stderr)
            # continua with this ion as the recombining/scattering ion
            if stage >= 2:
                try:
                    cont = ch.continuum(ions, float(ts[0]) if ts.size == 1 else ts, abundance=1.0)
                    cont.freeFree(wmid)
                    ff = np.atleast_2d(cont.FreeFree['intensity'])
                    try:
                        cont.freeBound(wmid)
                        fb = np.atleast_2d(cont.FreeBound['intensity'])
                    except Exception:
                        fb = np.zeros_like(ff)
                    for j, ij in enumerate(np.where(sel)[0]):
                        em[ij] += 4 * np.pi * (ff[j] + (fb[j] if fb.shape[0] > j else 0.0)) * dwl
                except Exception as e:
                    print('  %s continuum: %s' % (ions, str(e)[:80]), file=sys.stderr)
        tab[z] = em
        print('%-3s done in %.0f s; emissivity per n_e n_el at log T = %s: %s' % (
            el, time.time() - tstart, logt[[0, len(logt) // 2, -1]],
            ' '.join('%.3e' % s for s in em[[0, len(logt) // 2, -1]].sum(axis=1))), file=sys.stderr)
    with open(a.out, 'w') as f:
        f.write('# SuperLite CIE emissivities: CHIANTI database (Dere et al. 1997; Del Zanna et al. 2021)\n')
        f.write('# via ChiantiPy; per element and temperature, 4 pi eps per n_e n_el in each wavelength\n')
        f.write('# bin [erg cm^3 s^-1] (lines, free-free, free-bound, two-photon), CIE ion fractions\n')
        f.write('nel %d\nelements %s\n' % (len(tab), ' '.join(str(z) for z in tab)))
        f.write('nt %d\nlogt %s\n' % (len(logt), ' '.join('%.4f' % x for x in logt)))
        f.write('nw %d\nedges %s\n' % (nw, ' '.join('%.6e' % x for x in edges)))
        f.write('ne %.4e\n' % a.ne)
        for z, em in tab.items():
            for i, lt in enumerate(logt):
                f.write('%d %.4f\n' % (z, lt))
                f.write(' '.join('%.4e' % x for x in em[i]) + '\n')
    print('wrote', a.out, file=sys.stderr)


if __name__ == '__main__':
    main()
