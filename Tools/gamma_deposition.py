#!/usr/bin/env python3
"""Energy deposition of the 56Ni -> 56Co -> 56Fe decay in homologous ejecta.

Two methods:
 - gray absorption (default): kappa = 0.06 Y_e cm^2/g (Swartz, Sutherland &
   Harkness 1995), integrated along rays in the sphere,
       D(r) = kappa rho(r) (1/2) Int_-1^1 dmu Int_0^smax ds eps_g(r') exp(-tau(s));
   good to ~10 % once the ejecta are thin to gamma-rays (tau < ~2), but at
   larger depth it carries too much energy outward (the outer layers of the
   StaNdaRT toy01 model at 19.5 d get 2-7 times CMFGEN's deposition);
 - Monte Carlo (--mc): the gamma-ray lines of the two decays (NNDC
   intensities, the 511 keV annihilation photons of the 56Co positrons),
   Compton scattering with the Klein-Nishina cross-section and angular
   distribution on all electrons (bound and free, n_e = rho Y_e/m_u), the
   electron's recoil deposited where it scatters; photoabsorption
   2.0e-23 (Z/26)^5 (E/100 keV)^-3 cm^2 per atom (the photoelectric cross-
   section of iron at 100 keV, Z^5 E^-3 around it), with the mean Z^5 per
   gram of each zone (zmix), and photons below e_cut deposited where they
   arrive; static (no Doppler shifts, instantaneous transport).  The
   deposition is tallied along the paths (the expected loss per unit
   length, n_e sigma_KN <1 - e'/e> + kappa_pe, times the packet energy), so
   that thin outer zones that few packets scatter in still get a smooth
   estimate.
Positrons deposit their kinetic energy where they are made.  Decay energies
per decay: 56Ni 1.7206 MeV of gamma-rays, 56Co 3.6072 MeV of which 0.1158 MeV
is positron kinetic energy (the StaNdaRT toy models' values; the line lists
are scaled to these totals).

As a module: deposition(vedge, rho, xni0, t, ...) or deposition_mc(...) ->
per zone [erg/s/cm^3].  As a script it compares with the deposition of a
StaNdaRT reference code:

usage: gamma_deposition.py <model .dat> <edep_<model>_<code>.txt> <t1> [t2 ...]
       [--kappa 0.03] [--nmu 64] [--mc [--pe] [--npkt 200000] [--nz 100]] [--ye 0.5]
"""
import argparse
import numpy as np

DAY = 86400.0
MEV = 1.602176634e-6
AMU = 1.66053907e-24
MSUN = 1.989e33
T_NI = 6.075 / np.log(2.0)    # mean lives [d]
T_CO = 77.236 / np.log(2.0)
E_NI = 1.7206                 # MeV per decay, gamma-rays
E_CO = 3.6072                 # MeV per decay, all
E_POS = 0.1158                # MeV per 56Co decay, positron kinetic energy


def decay_power(xni0, t):
    """gamma and positron power per gram [erg/s/g] at t [d] for the initial
    56Ni mass fraction xni0."""
    n0 = xni0 / (56.0 * AMU)
    nni = n0 * np.exp(-t / T_NI)
    nco = n0 * T_CO / (T_CO - T_NI) * (np.exp(-t / T_CO) - np.exp(-t / T_NI))
    rni = nni / (T_NI * DAY)
    rco = nco / (T_CO * DAY)
    eg = (rni * E_NI + rco * (E_CO - E_POS)) * MEV
    ep = rco * E_POS * MEV
    return eg, ep


def deposition(vedge, rho, xni0, t, kappa=0.03, nmu=64, rout=None, vout=None):
    """Deposition [erg/s/cm^3] at the zone centres of a homologous model:
    vedge (n+1) zone edges [cm/s], rho (n) density at t [g/cm^3], xni0 (n)
    initial 56Ni mass fraction; t [d].  Returns (D_gamma, D_positron,
    gamma emissivity, positron emissivity), all per zone.  vout: the
    velocities [cm/s] at which D is wanted (default the zone centres)."""
    tt = t * DAY
    re = vedge * tt
    rc = 0.5 * (re[1:] + re[:-1])
    eg, ep = decay_power(xni0, t)
    epsg = eg * rho
    epsp = ep * rho
    rmax = re[-1]
    rq = rc if vout is None else np.asarray(vout) * tt
    mu, wmu = np.polynomial.legendre.leggauss(nmu)
    dr = np.min(np.diff(re))
    out = np.zeros(len(rq))
    for i, r0 in enumerate(rq):
        # distance to the outer edge along each direction (backwards along the ray)
        smax = r0 * mu + np.sqrt(np.maximum(rmax**2 - r0**2 * (1.0 - mu**2), 0.0))
        ns = int(np.max(smax) / (0.5 * dr)) + 2
        x = (np.arange(ns) + 0.5) / ns            # midpoints in [0,1] of s/smax
        s = smax[:, None] * x[None, :]
        r = np.sqrt(np.maximum(r0**2 + s**2 - 2.0 * r0 * s * mu[:, None], 0.0))
        k = np.clip(np.searchsorted(re, r, side='right') - 1, 0, len(rc) - 1)
        ds = smax[:, None] / ns
        dtau = kappa * rho[k] * ds
        tau = np.cumsum(dtau, axis=1) - 0.5 * dtau
        integ = np.sum(epsg[k] * np.exp(-tau) * ds, axis=1)
        out[i] = 0.5 * np.sum(wmu * integ)
    kq = np.clip(np.searchsorted(re, rq, side='right') - 1, 0, len(rc) - 1)
    return kappa * rho[kq] * out, epsp[kq], epsg[kq], epsp[kq]


# gamma-ray lines [keV, photons per decay] (NNDC; lines above 0.5 %)
LINES_NI = [(158.38, 0.988), (269.50, 0.365), (480.44, 0.365), (749.95, 0.495),
            (811.85, 0.860), (1561.80, 0.140)]
LINES_CO = [(846.77, 0.9994), (977.37, 0.0142), (1037.84, 0.1405), (1175.10, 0.0225),
            (1238.29, 0.6646), (1360.21, 0.0428), (1771.36, 0.1541), (2015.18, 0.0302),
            (2034.79, 0.0777), (2598.50, 0.1696), (3009.65, 0.0104), (3201.96, 0.0313),
            (3253.50, 0.0762), (3272.99, 0.0185), (511.0, 2 * 0.1957)]
MEC2 = 510.99895       # keV
RE2 = 7.9407877e-26    # r_e^2 [cm^2]


def sigma_kn(e):
    """Klein-Nishina cross-section per electron [cm^2] at photon energy e [keV]."""
    x = np.asarray(e, float) / MEC2
    l = np.log1p(2 * x)
    return 2 * np.pi * RE2 * ((1 + x) / x**2 * (2 * (1 + x) / (1 + 2 * x) - l / x)
                              + l / (2 * x) - (1 + 3 * x) / (1 + 2 * x)**2)


def kn_loss_table(n=240, nmu=2001):
    """Mean fractional energy loss of a Compton scattering, <1 - e'/e> over
    the Klein-Nishina angular distribution, on a log grid of e [keV]."""
    eg = np.logspace(0, np.log10(2e4), n)
    mu = np.linspace(-1.0, 1.0, nmu)
    x = eg[:, None] / MEC2
    p = 1.0 / (1.0 + x * (1.0 - mu[None, :]))
    w = p**2 * (p + 1.0 / p - (1.0 - mu[None, :]**2))     # dsigma/dmu up to a constant
    f = np.trapezoid(w * (1.0 - p), mu, axis=1) / np.trapezoid(w, mu, axis=1)
    return eg, f


def sample_kn(e, rng):
    """Energy ratio e'/e of a Compton scattering (Kahn's method), and the
    cosine of the scattering angle, for photon energies e [keV]."""
    x = e / MEC2
    n = len(x)
    out = np.empty(n)
    todo = np.arange(n)
    while todo.size:
        k = x[todo]
        r1, r2, r3 = rng.random((3, todo.size))
        a = (1 + 2 * k) / (9 + 2 * k)
        lo = r1 <= a
        eta = np.where(lo, 1 + 2 * k * r2, (1 + 2 * k) / (1 + 2 * k * r2))
        mu = 1 - (eta - 1) / k
        acc = np.where(lo, r3 <= 4 * (1 / eta - 1 / eta**2), r3 <= 0.5 * (mu**2 + 1 / eta))
        out[todo[acc]] = eta[acc]
        todo = todo[~acc]
    ratio = 1 / out
    return ratio, 1 - (out - 1) / x


def z5_per_gram(x):
    """sum_el X_el Z^5/(A_el m_u) [g^-1] for mass fractions x = {Z: X}."""
    a = {6: 12.011, 8: 15.999, 14: 28.085, 16: 32.06, 20: 40.078, 22: 47.867, 26: 55.845, 27: 58.933, 28: 58.693}
    return sum(xz * z**5 / (a.get(z, 2.0 * z) * AMU) for z, xz in x.items())


def deposition_mc(vedge, rho, xni0, t, ye=0.5, npkt=200000, ecut=30.0, seed=1, zmix=None):
    """Monte Carlo gamma-ray deposition [erg/s/cm^3] per zone of a homologous
    model (as deposition()): returns (D_gamma, D_positron).  zmix: per zone
    sum_el X_el Z^5/(A_el m_u) [g^-1] for the photoabsorption (z5_per_gram;
    default none)."""
    rng = np.random.default_rng(seed)
    if zmix is None:
        zmix = np.zeros(len(rho))
    hole = vedge[0] > 0
    if hole:
        # an empty zone inside the inner edge, crossed freely
        vedge = np.concatenate([[0.0], vedge])
        rho = np.concatenate([[0.0], rho])
        xni0 = np.concatenate([[0.0], xni0])
        zmix = np.concatenate([[0.0], zmix])
    kpe0 = 2.0e-23 / 26.0**5 * rho * zmix          # photoabsorption at 100 keV [1/cm]
    tt = t * DAY
    re = vedge * tt
    nz = len(rho)
    vol = 4.0 / 3.0 * np.pi * (re[1:]**3 - re[:-1]**3)
    eg, ep = decay_power(xni0, t)
    n0 = xni0 / (56.0 * AMU)
    rni = n0 * np.exp(-t / T_NI) / (T_NI * DAY) * rho * vol             # decays/s per zone
    rco = n0 * T_CO / (T_CO - T_NI) * (np.exp(-t / T_CO) - np.exp(-t / T_NI)) / (T_CO * DAY) * rho * vol
    eni = np.array([e for e, f in LINES_NI])
    fni = np.array([e * f for e, f in LINES_NI])
    fni = fni / fni.sum() * E_NI * 1e3                                     # keV per decay in each line
    eco = np.array([e for e, f in LINES_CO])
    fco = np.array([e * f for e, f in LINES_CO])
    fco = fco / fco.sum() * (E_CO - E_POS) * 1e3
# emission per (zone, line) [keV/s], packets of equal energy
    w = np.concatenate([np.outer(rni, fni), np.outer(rco, fco)], axis=1)
    eline = np.concatenate([eni, eco])
    wtot = w.sum()
    if wtot <= 0:
        return np.zeros(nz - int(hole)), (ep * rho)[int(hole):]
    cw = np.cumsum(w.ravel()) / wtot
    pick = np.minimum(np.searchsorted(cw, rng.random(npkt)), w.size - 1)
    iz = pick // len(eline)
    e = eline[pick % len(eline)].astype(float)
    wt = np.full(npkt, wtot / npkt)                                        # keV/s per packet
    r = np.cbrt(re[iz]**3 + rng.random(npkt) * (re[iz + 1]**3 - re[iz]**3))
    mu = 1 - 2 * rng.random(npkt)
    dep = np.zeros(nz)
    ne = rho * ye / AMU
    lte, ltf = kn_loss_table()
    alive = np.ones(npkt, bool)
    nloop = 0
    while alive.any():
        nloop += 1
        idx = np.nonzero(alive)[0]
        if nloop > 100000:
            raise RuntimeError('deposition_mc: packets stuck')
        k = np.clip(np.searchsorted(re, r[idx], side='right') - 1, 0, nz - 1)
        # distance to the shell boundaries
        rr, mm = r[idx], mu[idx]
        b2 = rr**2 * (1 - mm**2)
        dout = -rr * mm + np.sqrt(np.maximum(re[k + 1]**2 - b2, 0))
        din = np.full(idx.size, np.inf)
        hit = (mm < 0) & (re[k]**2 > b2)
        din[hit] = -rr[hit] * mm[hit] - np.sqrt(re[k[hit]]**2 - b2[hit])
        dbnd = np.minimum(dout, din)
        kcs = ne[k] * sigma_kn(e[idx])
        kpe = kpe0[k] * (e[idx] / 100.0)**-3
        kap = kcs + kpe
        dcol = -np.log(1 - rng.random(idx.size)) / np.maximum(kap, 1e-300)
        coll = dcol < dbnd
        # (across a boundary by a step resolved at the radius: a packet on the
        # boundary has dbnd = 0)
        d = np.where(coll, dcol, dbnd + 1e-9 * (dbnd + re[k + 1]))
        # new position and direction cosine after moving d
        rn = np.sqrt(rr**2 + d**2 + 2 * rr * d * mm)
        mun = (rr * mm + d) / np.maximum(rn, 1e-300)
        # expected deposition along the step, in the zone it starts in
        floss = np.interp(np.log(e[idx]), np.log(lte), ltf)
        np.add.at(dep, k, wt[idx] * np.where(coll, dcol, dbnd) * (kcs * floss + kpe))
        r[idx] = rn
        mu[idx] = np.clip(mun, -1, 1)
        esc = ~coll & (rn >= re[-1] * (1 - 1e-12))
        alive[idx[esc]] = False
        # photoabsorption: all of the packet's energy where it is absorbed
        pab = coll & (rng.random(idx.size) * kap < kpe)
        if pab.any():
            alive[idx[pab]] = False
        coll = coll & ~pab
        c = idx[coll]
        if c.size:
            kc = np.clip(np.searchsorted(re, r[c], side='right') - 1, 0, nz - 1)
            ratio, cth = sample_kn(e[c], rng)
            wt[c] *= ratio
            e[c] *= ratio
            # new direction: rotate by the scattering angle about the old one
            phi = 2 * np.pi * rng.random(c.size)
            m0 = mu[c]
            s0 = np.sqrt(np.maximum(1 - m0**2, 0))
            sth = np.sqrt(np.maximum(1 - cth**2, 0))
            mu[c] = np.clip(m0 * cth + s0 * sth * np.cos(phi), -1, 1)
            low = e[c] < ecut
            if low.any():
                np.add.at(dep, kc[low], wt[c[low]])
                alive[c[low]] = False
    dg = dep * 1e3 * (MEV / 1e6) / vol            # keV/s -> erg/s, per volume
    if hole:
        return dg[1:], (ep * rho)[1:]
    return dg, ep * rho


def rebin(vedge, dm, xni0, nz, zmix=None):
    """Zones merged to nz zones uniform in velocity (mass, 56Ni and the
    photoabsorption weight zmix conserved; zones split by volume)."""
    ve = np.linspace(vedge[0], vedge[-1], nz + 1)
    m = np.zeros(nz)
    mni = np.zeros(nz)
    mz = np.zeros(nz)
    if zmix is None:
        zmix = np.zeros(len(dm))
    for j in range(len(dm)):
        lo = np.clip(ve[:-1], vedge[j], vedge[j + 1])
        hi = np.clip(ve[1:], vedge[j], vedge[j + 1])
        f = (hi**3 - lo**3) / (vedge[j + 1]**3 - vedge[j]**3)
        m += f * dm[j]
        mni += f * dm[j] * xni0[j]
        mz += f * dm[j] * zmix[j]
    mm = np.maximum(m, 1e-300)
    return ve, m, mni / mm, mz / mm


def read_model(fname):
    rows = [l.split() for l in open(fname) if l.strip() and not l.startswith('#')]
    a = np.array(rows, dtype=float)
    # columns of the StaNdaRT toy models: stable IGE, 56Ni(t=0), IME as Ca, S, Si, unburnt C+O
    zm = (a[:, 4] + a[:, 5]) * 26.5**5 / (56.0 * AMU) + a[:, 16] * 20**5 / (40.08 * AMU) \
        + a[:, 17] * 16**5 / (32.06 * AMU) + a[:, 18] * 14**5 / (28.09 * AMU) + a[:, 8] * 7**5 / (14.0 * AMU)
    return dict(vel=a[:, 1] * 1e5, dm=a[:, 2] * MSUN, xni0=a[:, 5], zmix=zm)


def read_edep(fname):
    times, rows = None, []
    for l in open(fname):
        if l.startswith('#TIMES'):
            times = np.array(l.split(':')[1].split(), dtype=float)
        elif not l.startswith('#') and l.strip():
            rows.append(l.split())
    a = np.array(rows, dtype=float)
    return times, a[:, 0], a[:, 1:]


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('model')
    p.add_argument('edep')
    p.add_argument('times', type=float, nargs='+')
    p.add_argument('--kappa', type=float, default=0.03)
    p.add_argument('--nmu', type=int, default=64)
    p.add_argument('--mc', action='store_true', help='Monte Carlo transport with Compton scattering')
    p.add_argument('--npkt', type=int, default=200000)
    p.add_argument('--ye', type=float, default=0.5)
    p.add_argument('--nz', type=int, default=100, help='zones of the Monte Carlo grid')
    p.add_argument('--pe', action='store_true', help='photoabsorption with the model composition')
    a = p.parse_args()
    m = read_model(a.model)
    dv = 2.0 * m['vel'][0]
    ve = np.concatenate([m['vel'] - 0.5 * dv, [m['vel'][-1] + 0.5 * dv]])
    te, vr, ed = read_edep(a.edep)
    for t in a.times:
        it = int(np.argmin(np.abs(te - t)))
        tt = te[it] * DAY
        vol = 4.0 / 3.0 * np.pi * ((ve[1:] * tt)**3 - (ve[:-1] * tt)**3)
        rho = m['dm'] / vol
        if a.mc:
            vb, mb, xb, zb = rebin(ve, m['dm'], m['xni0'], a.nz, m['zmix'] if a.pe else None)
            rb = mb / (4.0 / 3.0 * np.pi * ((vb[1:] * tt)**3 - (vb[:-1] * tt)**3))
            dgz, dpz = deposition_mc(vb, rb, xb, te[it], a.ye, a.npkt, zmix=zb)
            vc = 0.5 * (vb[1:] + vb[:-1])
            d = np.interp(vr * 1e5, vc, dgz + dpz)
        else:
            dg, dp, eg, epz = deposition(ve, rho, m['xni0'], te[it], a.kappa, a.nmu, vout=vr * 1e5)
            d = dg + dp
        # totals: the reference on its velocity points (shells between midpoints)
        vm = vr * 1e5
        vb = np.concatenate([[vm[0] - 0.5 * (vm[1] - vm[0])], 0.5 * (vm[1:] + vm[:-1]),
                             [vm[-1] + 0.5 * (vm[-1] - vm[-2])]])
        sv = 4.0 / 3.0 * np.pi * ((vb[1:] * tt)**3 - (vb[:-1] * tt)**3)
        eg_all, ep_all = decay_power(m['xni0'], te[it])
        ptot = np.sum((eg_all + ep_all) * m['dm'])
        print('t = %7.2f d: decay power %.4e erg/s; deposited: here %.4e (%.3f), reference %.4e (%.3f)'
              % (te[it], ptot, np.sum(d * sv), np.sum(d * sv) / ptot, np.sum(ed[:, it] * sv),
                 np.sum(ed[:, it] * sv) / ptot))
        sel = [int(np.argmin(np.abs(vr - v))) for v in (2000, 5000, 8000, 11000, 14000, 18000, 22000, 26000)]
        print('   v [km/s]   ' + ' '.join('%9.0f' % vr[j] for j in sel))
        print('   here/ref   ' + ' '.join('%9.3f' % (d[j] / max(ed[j, it], 1e-300)) for j in sel))


if __name__ == '__main__':
    main()
