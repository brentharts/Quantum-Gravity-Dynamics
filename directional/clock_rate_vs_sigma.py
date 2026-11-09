#!/usr/bin/env python3
"""
clock_rate_vs_sigma.py -- the toy "clock rate = 1 / density" rule against the clock rate of the QGD Schwarzschild recipe.

The toy model GravitationalTimelineNetwork (spectre-monotile-py/research/toymodel_v1) gives every node a density rho >= 1 and
defines the local clock rate as 1/rho; a "well" adds  w / (1 + d^2)  to the density at distance d from its centre.  In QGD the
Schwarzschild recipe has sigma_t = sqrt(2GM/(c^2 r)) = sqrt(r_s/r) and g_tt = -(1 - sigma_t^2) (key-equations.md), so a static
observer's clock runs at  sqrt(1 - sigma_t^2) = sqrt(1 - r_s/r)  relative to infinity.

Two assumptions are needed to compare them, and both are stated here and not hidden: (a) the distance d of the toy well is read as
the radius r in units of r_s (in the toy it is a distance along a timeline, which is a different thing); (b) the density is a
free input of the toy, so any match is a statement about which density one would have to put in.

Results
-------
Exact 1.  The density that reproduces the QGD clock exactly is  rho_req(r) = (1 - r_s/r)^(-1/2).  It diverges at r = r_s.
          This is a definition of rho, not a derivation of the clock rate from the toy's dynamics.
Exact 2.  With x = GM/(c^2 r), the toy rule with rho = 1 + x gives 1/rho = 1 - x + x^2 - ..., and QGD / GR give
          sqrt(1 - 2x) = 1 - x - x^2/2 - ...: they agree to first order and differ at second order, by 3x^2/2.
Exact 3.  The Lorentzian well has rho - 1 = w/(1 + d^2) ~ w/d^2 at large d (exponent -2), whereas rho_req - 1 ~ r_s/(2r) (exponent -1).
          So the Lorentzian profile is not the density that the 1/r potential requires.
Exact 4.  The toy density is bounded by 1 + w, so its clock rate never drops below 1/(1 + w) > 0: the toy rule has no horizon for any
          finite well, whereas sqrt(1 - r_s/r) vanishes at r_s.
Numeric.  The maximum difference between the toy clock (with the first-order density 1 + r_s/(2r)) and the QGD clock over r/r_s in [1.5, 50].

    python3 directional/clock_rate_vs_sigma.py
"""
import numpy as np
import sympy as sp

x, w, d, rs, r = sp.symbols("x w d r_s r", positive=True)


def clock_qgd(r_over_rs):
    """Static-observer clock rate sqrt(1 - sigma_t^2) with sigma_t^2 = r_s / r."""
    return np.sqrt(1.0 - 1.0 / np.asarray(r_over_rs, dtype=float))


def density_required(r_over_rs):
    """rho_req = 1 / sqrt(1 - r_s/r)."""
    return 1.0 / clock_qgd(r_over_rs)


def clock_toy(rho):
    return 1.0 / np.asarray(rho, dtype=float)


def series_difference():
    """Exact: series of the two clocks in x = GM/(c^2 r)."""
    gr = sp.series(sp.sqrt(1 - 2 * x), x, 0, 4).removeO()
    toy = sp.series(1 / (1 + x), x, 0, 4).removeO()
    return sp.expand(gr), sp.expand(toy), sp.expand(toy - gr)


def tail_exponent(f, var, large=sp.oo):
    """Exponent n with f ~ C var^n as var -> infinity."""
    n = sp.limit(sp.log(f) / sp.log(var), var, large)
    return sp.simplify(n)


def main():
    gr, toy, diff = series_difference()
    print("clock rate series in x = GM/(c^2 r):")
    print("  QGD / GR  sqrt(1 - 2x) =", gr)
    print("  toy       1/(1 + x)    =", toy)
    print("  difference (toy - QGD) =", diff)
    print("\ntail exponents of rho - 1 (rho_req in r/r_s = u):")
    u = sp.symbols("u", positive=True)
    print("  Lorentzian well w/(1+d^2):", tail_exponent(w / (1 + d ** 2), d))
    print("  required (1 - 1/u)^(-1/2) - 1:", tail_exponent(sp.series((1 - 1 / u) ** sp.Rational(-1, 2) - 1, u, sp.oo, 3).removeO(), u))
    rr = np.linspace(1.5, 50, 500)
    rho1 = 1.0 + 1.0 / (2 * rr)          # first-order density  1 + r_s/(2r)
    print("\nmax |toy clock - QGD clock| with rho = 1 + r_s/(2r), r/r_s in [1.5, 50]: %.4f" % np.max(np.abs(clock_toy(rho1) - clock_qgd(rr))))
    print("toy clock rate floor for the default well (w = 4):", 1 / (1 + 4.0))


if __name__ == "__main__":
    main()
