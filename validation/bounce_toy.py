#!/usr/bin/env python3
"""
bounce_toy.py -- conditional bounce model for the l_Q^2 term (exploratory, not a claim of the theory).
=====================================================================================================

The repository has no covariant FRW form of sigma, so a Friedmann equation cannot be derived from the master metric.
The static term A/r^2 implies rho = -A/(8 pi G r^4).  ASSUMING it acts like a fluid  rho_q = -C a^-n  (the step r -> a is
an assumption) against matter of equation of state w, in units where 8 pi G / 3 = 1 and the matter density is 1 at a = 1:

    H^2 = a^(-3(1+w)) - C a^(-n)                                        (Friedmann)
    addot/a = -(1/2) [ (1+3w) a^(-3(1+w)) - (n-2) C a^(-n) ]            (rho_q conserved: w_q = n/3 - 1)

H = 0 exists iff n > 3(1+w), at a_b^(n-3(1+w)) = C, and addot/a = (n - 3(1+w)) C a_b^-n / 2 > 0 there automatically.
n = 4: dust bounces (a_b = C), radiation and stiffer matter do not.  A term  s*kappa*H^4  in  H^2 + s kappa H^4 = x  has
no bounce for either sign (s = +1: H^2 > 0 for all x > 0; s = -1: no real solution for x > 1/(4 kappa)).
"""
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp


def bounce_scale_factor(w, n, C):
    """Scale factor where H = 0, or None when n <= 3(1+w) (no bounce)."""
    p = 3.0 * (1.0 + w)
    if n <= p + 1e-12:
        return None
    return C ** (1.0 / (n - p))


def acceleration_sign_at_bounce(w, n, C=0.2):
    """Sign of addot at the bounce (+1 means it turns around); None when there is no bounce."""
    ab = bounce_scale_factor(w, n, C)
    if ab is None:
        return None
    X = ab ** (-n) * C
    return float(np.sign(-0.5 * ((1 + 3 * w) * X - (n - 2) * X)))


def integrate_bounce(w, n, C, tmax=12.0, a0=1.0, a_floor=1e-3):
    """Integrate addot = -(a/2)[(1+3w) a^-3(1+w) - (n-2) C a^-n] from a0, contracting (H < 0), until tmax or a < a_floor."""
    p = 3.0 * (1.0 + w)

    def rhs(t, y):
        a, v = y
        return [v, -0.5 * a * ((1 + 3 * w) * a ** (-p) - (n - 2) * C * a ** (-n))]

    def floor(t, y):
        return y[0] - a_floor
    floor.terminal, floor.direction = True, -1
    v0 = -a0 * np.sqrt(a0 ** (-p) - C * a0 ** (-n))
    sol = solve_ivp(rhs, (0, tmax), [a0, v0], events=floor, rtol=1e-11, atol=1e-13, max_step=1e-3)
    return sol.t, sol.y[0]


def h4_term_has_bounce():
    """Solve H^2 + s k H^4 = x for H^2.  Returns {s: [root_plus, root_minus]} for s = +1, -1."""
    k, x = sp.symbols("k x", positive=True)
    out = {}
    for s in (+1, -1):
        d = sp.sqrt(1 + 4 * s * k * x)
        out[s] = [(-1 + d) / (2 * s * k), (-1 - d) / (2 * s * k)]
    return out


if __name__ == "__main__":
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 4))
    for w, n, lab in [(0.0, 4.0, "dust, n=4"), (1 / 3, 4.0, "radiation, n=4"), (1 / 3, 6.0, "radiation, n=6")]:
        t, a = integrate_bounce(w, n, 0.2)
        ax.plot(t, a, label=lab)
        print(f"{lab:16s} a_b = {bounce_scale_factor(w, n, 0.2)}  min a = {a.min():.4g}")
    ax.set_xlabel("t"); ax.set_ylabel("a(t)"); ax.set_yscale("log"); ax.legend()
    fig.tight_layout(); fig.savefig(__file__.replace(".py", ".png"), dpi=120)
