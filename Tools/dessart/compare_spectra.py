#!/usr/bin/env python3
"""SuperLite formal-integral spectra (output.fi_spectrum: lambda [A], L_lambda)
against published CMFGEN spectra of Dessart et al. (Zenodo .fl files: lambda
[A], F_lambda at 1 kpc), in absolute L_lambda: both smoothed with the same
Gaussian (FWHM in A, as their figures), the band luminosities printed.
usage: compare_spectra.py out.png --pair label=superlite_spectrum,cmfgen.fl ...
                          [--fwhm 23.5] [--wl 1000,10000]"""
import argparse

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

KPC = 3.0857e21


def smooth(lam, f, fwhm):
    if fwhm <= 0:
        return f
    g = np.linspace(lam.min(), lam.max(), int((lam.max() - lam.min()) / (fwhm / 8)) + 1)
    fg = np.interp(g, lam, f)
    s = fwhm / 2.3548 / (g[1] - g[0])
    k = np.exp(-0.5 * (np.arange(-int(4 * s) - 1, int(4 * s) + 2) / s) ** 2)
    k /= k.sum()
    return np.interp(lam, g, np.convolve(fg, k, mode='same'))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('out')
    ap.add_argument('--pair', action='append', required=True)
    ap.add_argument('--fwhm', type=float, default=23.5)
    ap.add_argument('--wl', default='1000,10000')
    a = ap.parse_args()
    w1, w2 = [float(x) for x in a.wl.split(',')]
    bands = [(1000, 3000), (3000, 5000), (5000, 7000), (7000, 10000)]
    fig, ax = plt.subplots(len(a.pair), 1, figsize=(10, 3.6 * len(a.pair)), squeeze=False)
    for k, p in enumerate(a.pair):
        lab, rest = p.split('=', 1)
        fs, fc = rest.split(',')
        s = np.loadtxt(fs, comments='#')
        c = np.loadtxt(fc, comments='#')
        c = c[np.argsort(c[:, 0])]
        ls, Ls = s[:, 0], s[:, 1]
        lc, Lc = c[:, 0], 4 * np.pi * KPC ** 2 * c[:, 1]
        Ls, Lc = smooth(ls, Ls, a.fwhm), smooth(lc, Lc, a.fwhm)
        line = '%s:' % lab
        for b1, b2 in bands:
            ms, mc = (ls > b1) & (ls < b2), (lc > b1) & (lc < b2)
            if ms.sum() > 2 and mc.sum() > 2:
                es, ec = np.trapezoid(Ls[ms], ls[ms]), np.trapezoid(Lc[mc], lc[mc])
                line += '  L(%d-%d A) SL %.3e CMF %.3e (%.2f)' % (b1, b2, es, ec, es / ec)
        print(line)
        x = ax[k, 0]
        x.plot(lc, Lc / 1e39, color='#e87ba4', lw=1, label='CMFGEN (Dessart+22)')
        x.plot(ls, Ls / 1e39, color='#2a78d6', lw=1, label='SuperLite')
        x.set_xlim(w1, w2)
        x.set_ylabel(r'$L_\lambda$ [10$^{39}$ erg s$^{-1}$ $\AA^{-1}$]')
        x.set_title(lab, fontsize=10)
        x.legend(fontsize=8)
    ax[-1, 0].set_xlabel(r'$\lambda$ [$\AA$]')
    fig.tight_layout()
    fig.savefig(a.out, dpi=110)


if __name__ == '__main__':
    main()
