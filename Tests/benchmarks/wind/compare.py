"""Compare the SuperLite spectrum around H-alpha with a Sobolev profile.

usage: python3 compare.py run_dir profile.dat [--bin KMS] [--plot out.png]
       [--cont U1 U2] [--range UMIN UMAX]

profile.dat: columns u [km/s, toward the observer], F/F_c (sobolev.py).
The model flux is normalized by its mean over cont-band |u| (default
2500-6000 km/s) and rebinned to --bin km/s; errors from the packet counts.
"""
import argparse
import os

import numpy as np

C = 299792.458


def air(wl_cm):
    s2 = (1e-4 / wl_cm) ** 2
    return wl_cm / (1 + 8.34254e-5 + 2.406147e-2 / (130 - s2) + 1.5998e-4 / (38.9 - s2))


LAM0 = air(1.0 / (109677.58 * (1 / 4 - 1 / 9))) * 1e8      # H-alpha as in nltehmod, A


def read(d):
    with open(os.path.join(d, 'output.flx_grid')) as fh:
        fh.readline()
        edges = np.array(fh.readline().split(), float) * 1e8
    last = lambda f: [np.array(l.split(), float) for l in open(os.path.join(d, f)) if l.strip()][-1]
    return edges, last('output.flx_luminos'), last('output.flx_lumnum')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('run'); ap.add_argument('profile')
    ap.add_argument('--bin', type=float, default=100.0)
    ap.add_argument('--cont', nargs=2, type=float, default=[2500.0, 6000.0])
    ap.add_argument('--range', nargs=2, type=float, default=[-2500.0, 2500.0])
    ap.add_argument('--plot')
    a = ap.parse_args()
    edges, lum, num = read(a.run)
    ue = C * (1 - edges / LAM0)                 # bin edges in u (decreasing with wavelength)
    uc = 0.5 * (ue[1:] + ue[:-1])
    dlam = np.diff(edges)
    fl = lum / dlam
    cm = (np.abs(uc) > a.cont[0]) & (np.abs(uc) < a.cont[1])
    fc = fl[cm].sum() / cm.sum()
    # coarse bins in u
    be = np.arange(a.range[0], a.range[1] + 0.5 * a.bin, a.bin)
    fm, em, fa, ucen = [], [], [], []
    prof = np.loadtxt(a.profile)
    for lo, hi in zip(be[:-1], be[1:]):
        s = (uc >= lo) & (uc < hi)
        if not s.any():
            continue
        n = num[s].sum()
        f = fl[s].mean() / fc
        fm.append(f); em.append(f / np.sqrt(max(n, 1)))
        fa.append(np.interp(uc[s], prof[:, 0], prof[:, 1]).mean())
        ucen.append(0.5 * (lo + hi))
    fm, em, fa, ucen = map(np.array, (fm, em, fa, ucen))
    chi = (fm - fa) / em
    print(f'{a.run}: {len(fm)} bins of {a.bin:g} km/s, chi2/dof = {np.sum(chi**2)/len(fm):.2f}, '
          f'max |dev| = {np.abs(chi).max():.1f} sigma at u = {ucen[np.argmax(np.abs(chi))]:.0f} km/s; '
          f'continuum from {cm.sum()} bins')
    if a.plot:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 1, figsize=(7, 6), sharex=True, gridspec_kw=dict(height_ratios=[3, 1]))
        ax[0].plot(prof[:, 0], prof[:, 1], 'k-', lw=1, label='Sobolev formal solution')
        ax[0].errorbar(ucen, fm, em, fmt='o', ms=3, color='C3', label='SuperLite')
        ax[0].set_ylabel('F / F_c'); ax[0].legend(); ax[0].set_title(os.path.basename(os.path.normpath(a.run)))
        ax[1].plot(ucen, chi, 'o', ms=3, color='C3'); ax[1].axhline(0, color='k', lw=0.5)
        ax[1].set_ylabel('(MC-S)/sigma'); ax[1].set_xlabel('velocity toward observer (km/s)')
        fig.tight_layout(); fig.savefig(a.plot, dpi=120)


if __name__ == '__main__':
    main()
