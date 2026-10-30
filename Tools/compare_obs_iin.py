#!/usr/bin/env python3
# This file is part of SuperLite. SuperLite is released under the terms of the GNU GPLv3, see COPYING.
"""SuperLite IIn spectra (formal integral, output.fi_spectrum) against observed
spectra (rest frame; wavelength, flux columns, '#' comments): each normalized
to its median in a continuum window, the whole optical and an Halpha panel,
with the Halpha equivalent width and FWHM, the Halpha/Hbeta flux ratio (each
above a straight continuum) and the slope F(4400-4600)/F(6000-6200) printed.

usage: compare_obs_iin.py out.png --model label=path/output.fi_spectrum ...
                          --obs label=path/spectrum.dat ... [--norm 6000,6200]
                          [--res FWHM_kms]
--res convolves the model spectra with a Gaussian of that FWHM [km/s], the
resolution of the observed spectra (all measures are of the smoothed spectra).
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


HA = (6450, 6680, (6200, 6300), (6800, 6900))
HB = (4800, 4930, (4700, 4760), (4990, 5050))


def line(w, f, a, b, c1, c2):
    """equivalent width [A] (emission > 0) and flux above a straight
    continuum through the medians of the windows c1, c2; nan without them"""
    m1 = (w > c1[0]) & (w < c1[1]); m2 = (w > c2[0]) & (w < c2[1])
    if m1.sum() < 2 or m2.sum() < 2:
        return np.nan, np.nan
    x1, y1 = w[m1].mean(), np.median(f[m1]); x2, y2 = w[m2].mean(), np.median(f[m2])
    m = (w > a) & (w < b)
    c = y1 + (y2 - y1) * (w[m] - x1) / (x2 - x1)
    return np.trapezoid(f[m] / c - 1, w[m]), np.trapezoid(f[m] - c, w[m])


def ew(w, f):
    return line(w, f, *HA)[0]


def fwhm(w, f, l0=6562.8):
    """FWHM [km/s] of Halpha above the continuum of HA"""
    a, b, c1, c2 = HA
    m1 = (w > c1[0]) & (w < c1[1]); m2 = (w > c2[0]) & (w < c2[1])
    if m1.sum() < 2 or m2.sum() < 2:
        return np.nan
    x1, y1 = w[m1].mean(), np.median(f[m1]); x2, y2 = w[m2].mean(), np.median(f[m2])
    m = (w > 6300) & (w < 6820)
    y = f[m] - (y1 + (y2 - y1) * (w[m] - x1) / (x2 - x1))
    v = (w[m] / l0 - 1) * 2.99792458e5
    i = lo = hi = int(np.argmax(y))
    while lo > 0 and y[lo] > y[i] / 2:
        lo -= 1
    while hi < len(y) - 1 and y[hi] > y[i] / 2:
        hi += 1
    return v[hi] - v[lo]


def measures(w, f):
    """Halpha EW [A], FWHM [km/s], Halpha/Hbeta, F(4400-4600)/F(6000-6200)"""
    fa = line(w, f, *HA)[1]; fb = line(w, f, *HB)[1]
    mb = (w > 4400) & (w < 4600); mr = (w > 6000) & (w < 6200)
    slope = np.median(f[mb]) / np.median(f[mr]) if mb.sum() > 2 and mr.sum() > 2 else np.nan
    return ew(w, f), fwhm(w, f), fa / fb if fb > 0 else np.nan, slope


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
            w, f = smooth(w, f, a.res)
        m = (w > n1) & (w < n2)
        f = f / np.median(f[m])
        e, fw, dec, sl = measures(w, f)
        print(f'{lab}: Halpha EW {e:.0f} A, FWHM {fw:.0f} km/s, Halpha/Hbeta {dec:.2f}, '
              f'F(4500)/F(6100) {sl:.2f}')
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
