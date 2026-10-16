"""Sobolev line profile of a pure-scattering line in a spherical flow with
Q(mu) = mu^2 dv/dr + (1-mu^2) v/r > 0 everywhere (no resonance coupling between
different points), around an absorbing core of uniform intensity I_c.

  tau(r,mu) = tau0 (v/r)_cell / |Q(mu)|                 (as in the code)
  S(r) = beta_c I_c / beta,  beta = int_0^1 f dmu,  beta_c = 1/2 int_muc^1 f dmu,
  f(mu) = (1 - exp(-tau(mu)))/tau(mu),  muc = sqrt(1 - (RC/r)^2)
  I(p, u) = I_c e^-tau + S (1 - e^-tau)  (resonance in front of the core)
          = S (1 - e^-tau)              (p > RC),   I_c if occulted
  F(u)/F_c = int I 2 p dp / RC^2

The velocity field is the piecewise-linear one of input.str, as in SuperLite.
For a homologous flow this reduces to S = W(r) I_c (Jeffery & Branch 1990).

usage: python3 sobolev.py input.str tau0 out.dat [umax_kms] [nu]
"""
import sys

import numpy as np


def read_str(path):
    d = np.loadtxt(path, comments='#')
    r = np.r_[d[:, 0], d[-1, 1]]
    v = np.r_[d[:, 2], d[-1, 3]]
    return r, v


def profile(r, v, tau0, u, nmu=2001, nsub=41, nz=6000):
    n = len(r) - 1
    a = np.diff(v) / np.diff(r)                      # dv/dr per cell
    qv = (v[1:] + v[:-1]) / (r[1:] + r[:-1])         # (v/r)_cell
    rc, rout = r[0], r[-1]
    mu = np.linspace(0.0, 1.0, nmu)

    def tau_of(j, rr, m):
        vr = v[j] + a[j] * (rr - r[j])
        b = vr / rr
        q = np.abs(m ** 2 * a[j] + (1 - m ** 2) * b)
        q = np.maximum(q, 1e-3 * np.maximum(np.abs(a[j]), np.abs(b)))
        return tau0 * qv[j] / q

    # source function on nsub points per cell
    rs = np.empty((n, nsub)); ss = np.empty((n, nsub))
    for j in range(n):
        rs[j] = np.linspace(r[j], r[j + 1], nsub)
        for k, rr in enumerate(rs[j]):
            t = tau_of(j, rr, mu)
            f = np.where(t > 1e-8, -np.expm1(-t) / np.maximum(t, 1e-300), 1.0)
            beta = np.trapezoid(f, mu)
            muc = np.sqrt(max(0.0, 1 - (rc / rr) ** 2))
            m = mu >= muc
            betac = 0.5 * np.trapezoid(f[m], mu[m]) if m.sum() > 1 else 0.0
            ss[j, k] = betac / beta

    def s_at(rr):
        j = np.clip(np.searchsorted(r, rr, side='right') - 1, 0, n - 1)
        out = np.empty_like(rr)
        for jj in np.unique(j):
            sel = j == jj
            out[sel] = np.interp(rr[sel], rs[jj], ss[jj])
        return out, j

    p = np.r_[np.linspace(0, rc, 1200, endpoint=False), rc * (rout / rc) ** np.linspace(0, 1, 2400)][1:-1]
    flux = np.zeros_like(u)
    integ = np.zeros((len(p), len(u)))
    for ip, pp in enumerate(p):
        zmax = np.sqrt(rout ** 2 - pp ** 2)
        zmin = np.sqrt(rc ** 2 - pp ** 2) if pp < rc else -zmax
        z = np.linspace(zmin, zmax, nz)
        rr = np.hypot(pp, z)
        j = np.clip(np.searchsorted(r, rr, side='right') - 1, 0, n - 1)
        g = z / rr * (v[j] + a[j] * (rr - r[j]))      # projected velocity toward the observer
        ok = (u > g[0]) & (u < g[-1])
        I = np.where(pp < rc, 1.0, 0.0) * np.ones_like(u)
        if ok.any():
            zr = np.interp(u[ok], g, z)
            rr_ = np.hypot(pp, zr)
            S, jr = s_at(rr_)
            mres = zr / rr_
            t = tau_of(jr, rr_, mres)
            e = np.exp(-t)
            I[ok] = (e if pp < rc else 0.0) + S * (1 - e)
        integ[ip] = I * 2 * pp
    flux = np.trapezoid(integ, p, axis=0) / rc ** 2
    return flux


if __name__ == '__main__':
    path, tau0, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    umax = float(sys.argv[4]) if len(sys.argv) > 4 else 2500.0
    nu = int(sys.argv[5]) if len(sys.argv) > 5 else 251
    r, v = read_str(path)
    u = np.linspace(-umax, umax, nu) * 1e5
    f = profile(r, v, tau0, u)
    np.savetxt(out, np.c_[u / 1e5, f], header='u = velocity toward observer [km/s], F/F_c')
    print(f'min {f.min():.4f} at {u[np.argmin(f)]/1e5:.0f} km/s, max {f.max():.4f} at {u[np.argmax(f)]/1e5:.0f} km/s')
