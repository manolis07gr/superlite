"""Check sobolev_nl.py against Rybicki & Hummer (1978, ApJ 219, 654), Table 2:
a pure-scattering line (Ic = 1, B = eps = 0) in the decelerating winds
v = v0 (rc/r)^alpha, alpha = 1/2 and 1, out to 100 rc, tau0 = 100 (their
tau(mu) = tau0 (rc/r)^(1-2 alpha)/|1 - (1+alpha) mu^2|): the generalized
source function S and the disconnected one S0 (the core attenuated by the
resonances upstream, without their emission) at log r/rc = 0, 0.5, 1, 1.5, 2.
The velocity is linear and K constant in each of ncell logarithmic cells, as
in SuperLite.

usage: python3 rh78_table2.py [ncell (100)]
"""
import sys
import time

import numpy as np

import sobolev_nl as SN

# Table 2, block tau0 = 100: S_num and S0 (two significant digits)
S_RH = {0.5: [5.0e-1, 1.5e-2, 2.6e-3, 3.7e-4, 5.1e-5], 1.0: [5.2e-1, 4.8e-2, 9.3e-3, 1.7e-3, 3.1e-4]}
S0_RH = {0.5: [5.0e-1, 1.3e-2, 1.9e-3, 2.1e-4, 2.1e-5], 1.0: [5.0e-1, 3.5e-2, 4.0e-3, 4.1e-4, 4.1e-5]}


def main():
    ncell = int(sys.argv[1]) if len(sys.argv) > 1 else 100
    worst = 0.0
    for alpha in (0.5, 1.0):
        r = 100.0**(np.arange(ncell + 1)/ncell)
        v = r**-alpha
        fl = SN.Flow(r, v, 100.0)
        fl.k = 100.0*fl.qv*(0.5*(r[1:] + r[:-1]))**(2*alpha - 1)
        t = time.time()
        rs, S, bbar, sloc, jinc, s0 = SN.source_function(fl, nsub=4, nmu=400, nz=4000)
        lr = np.log10(rs)
        print('alpha %.1f (%.0f s)   log r    S        RH78     S0       RH78' % (alpha, time.time() - t))
        for i, x in enumerate([0.0, 0.5, 1.0, 1.5, 2.0]):
            xx = min(max(x, lr[0]), lr[-1])
            a = np.exp(np.interp(xx, lr, np.log(S)))
            b = np.exp(np.interp(xx, lr, np.log(s0)))
            worst = max(worst, abs(a/S_RH[alpha][i] - 1), abs(b/S0_RH[alpha][i] - 1))
            print('                      %.1f  %.3e  %.1e  %.3e  %.1e' % (x, a, S_RH[alpha][i], b, S0_RH[alpha][i]))
    print('largest relative difference %.3f' % worst)


if __name__ == '__main__':
    main()
