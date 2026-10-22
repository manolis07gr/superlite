"""Generalized Sobolev solution (Rybicki & Hummer 1978) of a pure-scattering
line in a spherical flow of any velocity law, around an absorbing core of
uniform intensity I_c = 1: the source function with the nonlocal coupling
between the points of a ray with the same projected velocity, and the
emergent profile.

  tau(r,mu) = K(r)/|Q(r,mu)|,  Q = mu^2 dv/dr + (1-mu^2) v/r,  K = tau0 (v/r)_cell
  beta(r,mu) = (1 - exp(-tau))/tau
  J(r) = (1 - beta_bar) S(r) + 1/2 int_-1^1 beta(r,mu) I_inc(r,mu) dmu
  S = J (pure scattering):  S(r) beta_bar(r) = 1/2 int beta I_inc dmu

I_inc(r,mu) is the intensity that reaches the resonance at r along the ray
of direction mu: from the core (I_c) or the outer boundary (0), through every
resonance upstream on the ray, I -> I e^-tau' + S(r') (1 - e^-tau'), at the
points where the projected velocity v(r') z'/r' equals v(r) mu.  In a flow
with Q > 0 everywhere there are none, and S = beta_c I_c/beta_bar (sobolev.py).
The velocity field and K are those of SuperLite's cells (piecewise-linear v,
K constant per cell, |Q| floored at 1e-3 max(|dv/dr|, v/r) as in the code).

usage: python3 sobolev_nl.py input.str tau0 out_prefix [umax_kms] [nu]
  writes out_prefix.S (r, S/I_c, beta_bar, S of the local Sobolev solution,
  the plain angle mean of the incident intensity, Rybicki & Hummer's S0)
  and out_prefix.dat (u [km/s] toward the observer, F/F_c)
"""
import sys

import numpy as np

from sobolev import read_str


class Flow:
    def __init__(self, r, v, tau0):
        self.r, self.v = r, v
        self.n = len(r) - 1
        self.a = np.diff(v)/np.diff(r)
        self.qv = (v[1:] + v[:-1])/(r[1:] + r[:-1])
        self.k = tau0*self.qv
        self.rc, self.rout = r[0], r[-1]

    def cell(self, rr):
        return np.clip(np.searchsorted(self.r, rr, side='right') - 1, 0, self.n - 1)

    def vel(self, rr):
        j = self.cell(rr)
        return self.v[j] + self.a[j]*(rr - self.r[j]), j

    def tau(self, rr, mu):
        vr, j = self.vel(rr)
        b = vr/rr
        q = np.abs(mu**2*self.a[j] + (1 - mu**2)*b)
        q = np.maximum(q, 1e-3*np.maximum(np.abs(self.a[j]), np.abs(b)))
        return self.k[j]/q


def beta_of(t):
    return np.where(t > 1e-8, -np.expm1(-t)/np.maximum(t, 1e-300), 1.0)


def crossings(f, z, w, w0, ztol):
    """Points on the ray (samples z, projected velocity w) where w = w0,
    by linear interpolation of the sign changes and two secant steps on
    the exact w(z) = f(z); those within ztol of the end of the samples are
    left out (the local resonance)."""
    d = w - w0
    s = np.nonzero(np.sign(d[:-1])*np.sign(d[1:]) < 0)[0]
    if s.size == 0:
        return np.empty(0)
    za, zb, da, db = z[s], z[s + 1], d[s], d[s + 1]
    zr = za - da*(zb - za)/(db - da)
    for _ in range(2):
        dr = f(zr) - w0
        left = np.sign(dr) == np.sign(da)
        za = np.where(left, zr, za)
        da = np.where(left, dr, da)
        zb = np.where(left, zb, zr)
        db = np.where(left, db, dr)
        zr = za - da*(zb - za)/np.where(db != da, db - da, 1.0)
    return zr[np.abs(zr - z[-1]) > ztol]


def ray_points(fl, p, z0, z1, w0, nz):
    """Resonances with projected velocity w0 on the ray of impact parameter p
    from z0 to z1 (increasing z), with their radius and direction cosine."""
    z = np.linspace(z0, z1, nz)

    def f(zz):
        rr = np.hypot(p, zz)
        return fl.vel(rr)[0]*zz/rr

    zr = crossings(f, z, f(z), w0, 1e-9*fl.rout)
    rr = np.hypot(p, zr)
    return rr, zr/rr


def source_function(fl, nsub=8, nmu=400, nz=3000, eps=0.0, bth=0.0):
    """S(r)/I_c on nsub points per cell (one linear system with the nonlocal
    coupling), and the local Sobolev solution.  With a thermalization
    probability eps and a thermal source bth (in units of I_c), as Rybicki &
    Hummer's two-level atom: S = (1-eps) J + eps B, J = (1-beta_bar) S + <beta I_inc>,
    S (eps + (1-eps) beta_bar) = (1-eps) <beta I_inc> + eps B."""
    rs = np.concatenate([np.linspace(fl.r[j], fl.r[j + 1], nsub, endpoint=(j == fl.n - 1))
                         for j in range(fl.n)])
    rs = np.clip(rs, fl.rc*(1 + 1e-9), fl.rout*(1 - 1e-9))
    xg, wg = np.polynomial.legendre.leggauss(nmu//2)
    xc, wc = np.polynomial.legendre.leggauss(max(nmu//4, 8))

    def mugrid(muc):
        # Gauss-Legendre on [-1,0], [0,mu_c] and the cone of the core [mu_c,1]
        # (far out it is narrow: 5e-5 wide at 100 R_c), 1/2 int dmu
        a = 0.5*(xc + 1)
        m = np.r_[-0.5*(xg + 1)[::-1], muc*a, muc + (1 - muc)*a]
        w = np.r_[0.5*wg[::-1], muc*0.5*wc, (1 - muc)*0.5*wc]/2
        return m, w
    # per (radius, mu): the upstream resonances (their S nodes and weights) and the
    # attenuated boundary intensity; I_inc = I0 + sum_j c_j S(r_j)
    rows, rowj = [], []
    bbar = np.empty(rs.size)
    sloc = np.empty(rs.size)
    for i, r0 in enumerate(rs):
        v0 = fl.vel(np.array([r0]))[0][0]
        muc = np.sqrt(max(0.0, 1 - (fl.rc/r0)**2))
        mu, wmu = mugrid(muc)
        t0 = fl.tau(np.full(mu.size, r0), mu)
        b0 = beta_of(t0)
        bbar[i] = np.sum(wmu*b0)
        p = r0*np.sqrt(1 - mu**2)
        z0 = r0*mu
        coef, coefj = [], []
        ibnd, jbnd = 0.0, 0.0
        for k in range(mu.size):
            if p[k] < fl.rc and z0[k] > 0:
                zs, i_start = np.sqrt(fl.rc**2 - p[k]**2), 1.0
            else:
                zs, i_start = -np.sqrt(fl.rout**2 - p[k]**2), 0.0
            if z0[k] - zs <= 0:
                continue
            rr, mm = ray_points(fl, p[k], zs, z0[k], v0*mu[k], nz)
            ii = i_start
            c, cj = {}, {}
            if rr.size:
                o = np.argsort(mm*rr)                             # increasing z
                rr, mm = rr[o], mm[o]
                t = fl.tau(rr, mm)
                e = np.exp(-t)
                # I_inc = i_start prod(e) + sum_j S_j (1-e_j) prod_{l>j} e_l
                tail = np.r_[np.cumprod(e[::-1])[::-1][1:], 1.0]
                ii = i_start*np.prod(e)
                for rj, wj in zip(rr, (1 - e)*tail):
                    # linear interpolation of S between nodes
                    m = np.searchsorted(rs, rj)
                    m = min(max(m, 1), rs.size - 1)
                    f = (rj - rs[m - 1])/(rs[m] - rs[m - 1])
                    c[m - 1] = c.get(m - 1, 0.0) + wmu[k]*b0[k]*wj*(1 - f)
                    c[m] = c.get(m, 0.0) + wmu[k]*b0[k]*wj*f
                    cj[m - 1] = cj.get(m - 1, 0.0) + wmu[k]*wj*(1 - f)
                    cj[m] = cj.get(m, 0.0) + wmu[k]*wj*f
            ibnd += wmu[k]*b0[k]*ii
            jbnd += wmu[k]*ii
            coef.append(c)
            coefj.append(cj)
        for co, bnd, out in ((coef, ibnd, rows), (coefj, jbnd, rowj)):
            cc = {}
            for c in co:
                for m, x in c.items():
                    cc[m] = cc.get(m, 0.0) + x
            out.append((bnd, cc))
        # the local solution: no nonlocal coupling
        sloc[i] = ((1 - eps)*np.sum((wmu*b0)[mu >= muc]) + eps*bth)/(eps + (1 - eps)*bbar[i])
    # S (eps + (1-eps) bbar) = (1-eps) (ibnd + sum c S) + eps B  (solve directly)
    n = rs.size
    A = np.zeros((n, n))
    b = np.empty(n)
    for i, (ibnd, cc) in enumerate(rows):
        A[i, i] += eps + (1 - eps)*bbar[i]
        for m, x in cc.items():
            A[i, m] -= (1 - eps)*x
        b[i] = (1 - eps)*ibnd + eps*bth
    S = np.linalg.solve(A, b)
    # the plain angle mean of the incident intensity, 1/2 int I_inc dmu
    jinc = np.array([bnd + sum(x*S[m] for m, x in cc.items()) for bnd, cc in rowj])
    # Rybicki & Hummer's disconnected approximation S0: the core seen through
    # the resonances upstream (their attenuation), without their emission
    s0 = b/np.diag(A)
    return rs, S, bbar, sloc, jinc, s0


def profile(fl, rs, S, u, npc=600, npo=1600, nz=4000):
    """F(u)/F_c: rays of impact parameter p, every resonance on each."""
    p = np.r_[np.linspace(0, fl.rc, npc, endpoint=False),
              fl.rc*(fl.rout/fl.rc)**np.linspace(0, 1, npo)][1:-1]
    integ = np.zeros((p.size, u.size))
    for ip, pp in enumerate(p):
        zmax = np.sqrt(fl.rout**2 - pp**2)
        zmin = np.sqrt(fl.rc**2 - pp**2) if pp < fl.rc else -zmax
        z = np.linspace(zmin, zmax, nz)
        rr = np.hypot(pp, z)
        w = fl.vel(rr)[0]*z/rr
        d = w[None, :] - u[:, None]
        sc = np.sign(d[:, :-1])*np.sign(d[:, 1:]) < 0
        I = np.full(u.size, 1.0 if pp < fl.rc else 0.0)
        iu, iz = np.nonzero(sc)
        if iu.size:
            za, zb = z[iz], z[iz + 1]
            da, db = d[iu, iz], d[iu, iz + 1]
            zr = za - da*(zb - za)/(db - da)
            r_ = np.hypot(pp, zr)
            t = fl.tau(r_, zr/r_)
            e = np.exp(-t)
            s_ = np.interp(r_, rs, S)
            # resonances in order of increasing z for each u (np.nonzero is row-major)
            for q in range(iu.size):
                I[iu[q]] = I[iu[q]]*e[q] + s_[q]*(1 - e[q])
        integ[ip] = I*2*pp
    return np.trapezoid(integ, p, axis=0)/fl.rc**2


if __name__ == '__main__':
    path, tau0, out = sys.argv[1], float(sys.argv[2]), sys.argv[3]
    umax = float(sys.argv[4]) if len(sys.argv) > 4 else 2500.0
    nu = int(sys.argv[5]) if len(sys.argv) > 5 else 251
    r, v = read_str(path)
    fl = Flow(r, v, tau0)
    rs, S, bbar, sloc, jinc, s0 = source_function(fl)
    np.savetxt(out + '.S', np.c_[rs, S, bbar, sloc, jinc, s0],
               header='r [cm], S/I_c (generalized Sobolev; = <beta I_inc>/beta_bar), beta_bar, '
                      'S/I_c (local Sobolev), <I_inc>/I_c (plain angle mean of the incident intensity), '
                      'S0/I_c (disconnected: the core attenuated by the upstream resonances, RH78)')
    u = np.linspace(-umax, umax, nu)*1e5
    f = profile(fl, rs, S, u)
    np.savetxt(out + '.dat', np.c_[u/1e5, f], header='u = velocity toward observer [km/s], F/F_c')
    print(f'S/I_c from {S.min():.4f} to {S.max():.4f}; local Sobolev {sloc.min():.4f} to {sloc.max():.4f}')
    print(f'profile min {f.min():.4f} at {u[np.argmin(f)]/1e5:.0f} km/s, max {f.max():.4f} at {u[np.argmax(f)]/1e5:.0f} km/s')
