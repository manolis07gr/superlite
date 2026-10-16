"""Analytic Sobolev P-Cygni profile (pure scattering, uniform tau, homologous,
sharp absorbing photosphere): S(r) = W(r) I_ph, formal integral over impact
parameters (Jeffery & Branch 1990).  Velocities in units of v_ph."""
import numpy as np

def profile(u, tau, vmax):
    """u: observed velocity toward the observer (blueshift > 0), in v_ph units.
    Returns F/F_continuum."""
    out = np.empty_like(u)
    ps = np.linspace(0, vmax, 4001)[1:]
    dp = ps[1] - ps[0]
    for i, z in enumerate(u):
        F = 0.0
        for p in ps - 0.5 * dp:
            r = np.hypot(p, z)
            inres = (r >= 1.0) and (r <= vmax)
            if p < 1.0:
                zph = np.sqrt(1 - p * p)        # photosphere surface along the ray
                if inres and z > zph:
                    W = 0.5 * (1 - np.sqrt(1 - 1 / r**2))
                    I = np.exp(-tau) + W * (1 - np.exp(-tau))
                else:
                    I = 1.0
            else:
                if inres:
                    W = 0.5 * (1 - np.sqrt(1 - 1 / r**2))
                    I = W * (1 - np.exp(-tau))
                else:
                    I = 0.0
            F += I * 2 * p * dp
        out[i] = F               # continuum: integral of 2 p dp over p<1 = 1
    return out

if __name__ == "__main__":
    u = np.linspace(-3.2, 3.2, 129)
    f = profile(u, 3.0, 3.0)
    np.savetxt("analytic_tau3.dat", np.c_[u, f], header="u=v_toward_observer/v_ph  F/Fc")
    print("min", f.min(), "at u=", u[np.argmin(f)], " max", f.max(), "at u=", u[np.argmax(f)])
