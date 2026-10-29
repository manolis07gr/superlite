#!/usr/bin/env python3
# This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
"""SuperLite IIn spectra (formal integral, output.fi_spectrum) against observed
spectra (rest frame; wavelength, flux columns, '#' comments): each normalized
to its median in a continuum window, the whole optical and an Halpha panel,
with the Halpha equivalent widths printed.

usage: compare_obs_iin.py out.png --model label=path/output.fi_spectrum ...
                          --obs label=path/spectrum.dat ... [--norm 6000,6200]
                          [--res FWHM_kms]
--res convolves the model spectra with a Gaussian of that FWHM [km/s], the
resolution of the observed spectra (the equivalent widths do not change).
"""
import argparse

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

COLORS = ['#2a78d6', '#1baf7a', '#6250d6', '#e87ba4']   # models
OBSCOL = ['#eb6834', '#eda100', '#e34948', '#008300']   # observations


def air(lam):
    """vacuum to air above 2000 A (SuperLite's spectra are in air there)"""
    s2 = (1e4 / lam)**2
    n = 1 + 0.0000834254 + 0.02406147 / (130 - s2) + 0.00015998 / (38.9 - s2)
    return np.where(lam > 2000, lam / n, lam)


def ew(w, f, a=6450, b=6680, c1=(6200, 6300), c2=(6800, 6900)):
    m1 = (w > c1[0]) & (w < c1[1]); m2 = (w > c2[0]) & (w < c2[1])
    if m1.sum() < 2 or m2.sum() < 2:
        return np.nan
    x1, y1 = w[m1].mean(), np.median(f[m1]); x2, y2 = w[m2].mean(), np.median(f[m2])
    m = (w > a) & (w < b)
    return np.trapezoid(f[m] / (y1 + (y2 - y1) * (w[m] - x1) / (x2 - x1)) - 1, w[m])


def smooth(w, f, fwhm):
    """Gaussian of FWHM [km/s] on a uniform grid in ln(wavelength)"""
    dl = fwhm / 2.99792458e5 / 5.
    lg = np.arange(np.log(w[0]), np.log(w[-1]), dl)
    fg = np.interp(lg, np.log(w), f)
    sig = fwhm / 2.99792458e5 / 2.3548 / dl
    k = np.arange(-int(4 * sig) - 1, int(4 * sig) + 2)
    ker = np.exp(-0.5 * (k / sig)**2)
    return np.exp(lg), np.convolve(fg, ker / ker.sum(), mode='same')


def load(spec, model):
    lab, path = spec.split('=', 1)
    d = np.loadtxt(path, comments='#')
    w, f = d[:, 0], d[:, 1]
    o = np.argsort(w)
    return lab, w[o], f[o]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('out')
    ap.add_argument('--model', action='append', default=[])
    ap.add_argument('--obs', action='append', default=[])
    ap.add_argument('--norm', default='6000,6200')
    ap.add_argument('--vacuum-obs', action='store_true', help='observed wavelengths in vacuum')
    ap.add_argument('--res', type=float, default=0., help='model smoothing FWHM [km/s]')
    a = ap.parse_args()
    n1, n2 = (float(x) for x in a.norm.split(','))
    fig, ax = plt.subplots(2, 1, figsize=(12, 8))
    sets = [(load(s, True), COLORS[i % 4], 1.3) for i, s in enumerate(a.model)] + \
           [(load(s, False), OBSCOL[i % 4], 1.0) for i, s in enumerate(a.obs)]
    for k, ((lab, w, f), col, lw) in enumerate(sets):
        if k >= len(a.model) and a.vacuum_obs:
            w = air(w)
        if k < len(a.model) and a.res > 0:
            e0 = ew(w, f)
            w, f = smooth(w, f, a.res)
        m = (w > n1) & (w < n2)
        f = f / np.median(f[m])
        e = ew(w, f) if not (k < len(a.model) and a.res > 0) else e0
        print(f'{lab}: Halpha EW {e:.1f} A')
        for j, (x1, x2) in enumerate([(3500, 9000), (6250, 6900)]):
            mm = (w > x1) & (w < x2)
            ax[j].plot(w[mm], f[mm], color=col, lw=lw, label=f'{lab} (Hα EW {e:.0f} Å)' if j == 0 else lab)
    for j, (x1, x2) in enumerate([(3500, 9000), (6250, 6900)]):
        ax[j].set_xlim(x1, x2)
        ax[j].set_ylabel(f'F$_\\lambda$ / median {n1:.0f}–{n2:.0f} Å')
        ax[j].grid(color='#e4e3df', lw=0.6)
        for s in ('top', 'right'):
            ax[j].spines[s].set_visible(False)
    ax[0].legend(frameon=False, fontsize=8)
    ax[1].set_xlabel('rest wavelength [Å] (air)')
    fig.tight_layout()
    fig.savefig(a.out, dpi=130)


if __name__ == '__main__':
    main()
