#!/usr/bin/env python3
"""
two_source_clocks.py -- the clock rate near two masses, three ways: additive density (the toy), linearised GR (Newtonian), and the QGD sigma-sum.

Single source: sigma_t^2 = r_s/r and the static clock is sqrt(1 - sigma_t^2).  For two sources key-equations.md (section 13) writes
g_tt = -[1 - (sigma_1 + sigma_2)^2], with sigma_a = sqrt(r_sa / r_a), and boxes the cross term 2 sigma_1 sigma_2 ~ r^-1.

  (i)   toy:    rho = 1 + x_1 + x_2, clock = 1/rho, x_a = G M_a / (c^2 r_a) = r_sa / (2 r_a).
  (ii)  linearised GR (static sources add in the Newtonian potential): clock = 1 - x_1 - x_2 + O(x^2).
  (iii) QGD sigma-sum: clock = sqrt(1 - (sigma_1 + sigma_2)^2) = 1 - x_1 - x_2 - 2 sqrt(x_1 x_2) + ...

Results
Exact 1.  (i) and (ii) agree to first order in G for any number of sources.  (iii) differs from them at FIRST order, by the cross term 2 sqrt(x_1 x_2),
          because sigma_a is proportional to sqrt(G M_a) and the cross term is proportional to sqrt(G M_1 G M_2) = G sqrt(M_1 M_2).
Exact 2.  Far from two masses (r_1 = r_2 = R >> separation) the sigma-sum clock is that of a single mass (sqrt M_1 + sqrt M_2)^2, against M_1 + M_2.
          For M_1 = M_2 this is twice the Newtonian mass.
Exact 3.  The cross term is not harmonic: Laplacian(1/sqrt(r_1 r_2)) != 0 off the sources, whereas Laplacian(1/r_a) = 0.  A weak-field static g_tt must have a harmonic
          potential in vacuum, so the sigma-sum cannot be a vacuum weak-field solution at first order in G.
Numeric.  (tests/test_directional_nbody.py)  For the sigma-sum metric the Ricci residual of tests/test_gaps.py halves when the masses halve (first order in G);
          for the linear sum of Kerr-Schild terms it quarters (second order), the scaling asserted in tests/test_gaps.py.  README.md says "order (GM)^2" for the sigma-sum.

    python3 directional/two_source_clocks.py
"""
import numpy as np
import sympy as sp

x1, x2, M1, M2, R = sp.symbols("x1 x2 M1 M2 R", positive=True)
X, Y, Z = sp.symbols("X Y Z", real=True)


def clock_toy(x_list):
    return 1.0 / (1.0 + sum(x_list))


def clock_linear_gr(x_list):
    return 1.0 - sum(x_list)


def clock_sigma_sum(x_list):
    """sigma_a^2 = 2 x_a, so sigma_total = sum_a sqrt(2 x_a), clock = sqrt(1 - sigma_total^2)."""
    s = sum(np.sqrt(2.0 * np.asarray(x)) for x in x_list)
    return np.sqrt(1.0 - s ** 2)


def first_order_cross_term():
    """Exact 1: coefficient of the cross term in the clock expansion, with x_a = eps^2 * (GM_a/r_a) scaling, in sympy."""
    eps = sp.symbols("epsilon", positive=True)
    s = sp.sqrt(2 * x1) + sp.sqrt(2 * x2)
    clock = sp.sqrt(1 - eps ** 2 * s ** 2)
    ser = sp.series(clock, eps, 0, 3).removeO()
    return sp.simplify(sp.expand(ser))


def far_field_mass_ratio():
    """Exact 2: effective mass of the sigma-sum over M1 + M2 when r1 = r2 = R."""
    sigma1_sq = M1 / R
    sigma2_sq = M2 / R
    total = (sp.sqrt(sigma1_sq) + sp.sqrt(sigma2_sq)) ** 2
    return sp.simplify(total / ((sigma1_sq + sigma2_sq)))


def laplacian_cross_term_at(point=(1.3, 0.4, 0.7), a=(-1.0, 0.0, 0.0), b=(1.0, 0.0, 0.0)):
    """Exact 3: Laplacian of 1/sqrt(r1 r2) at a point off both sources (sources at a and b)."""
    r1 = sp.sqrt((X - a[0]) ** 2 + (Y - a[1]) ** 2 + (Z - a[2]) ** 2)
    r2 = sp.sqrt((X - b[0]) ** 2 + (Y - b[1]) ** 2 + (Z - b[2]) ** 2)
    f = 1 / sp.sqrt(r1 * r2)
    lap = sum(sp.diff(f, v, 2) for v in (X, Y, Z))
    lap1 = sum(sp.diff(1 / r1, v, 2) for v in (X, Y, Z))
    sub = dict(zip((X, Y, Z), point))
    return float(lap.subs(sub)), float(sp.simplify(lap1.subs(sub)))


def main():
    print("clock series in eps (x_a -> eps^2 x_a):", first_order_cross_term())
    print("far-field mass ratio (sigma-sum / Newtonian), r1 = r2:", far_field_mass_ratio(), "-> equal masses:", far_field_mass_ratio().subs({M1: 1, M2: 1}))
    lap_cross, lap_single = laplacian_cross_term_at()
    print("Laplacian of 1/sqrt(r1 r2) off the sources: %.4f ; of 1/r1: %.2e" % (lap_cross, lap_single))
    for x in (0.001, 0.01):
        print("x1 = x2 = %.3f: toy %.6f, linear GR %.6f, sigma-sum %.6f" % (x, clock_toy([x, x]), clock_linear_gr([x, x]), clock_sigma_sum([x, x])))


if __name__ == "__main__":
    main()
