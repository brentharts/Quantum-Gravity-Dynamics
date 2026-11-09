"""Checks of directional/axis_timelines.py and directional/clock_rate_vs_sigma.py.

The algebra for the QGD side is redone in sympy (Kerr-Schild form, as in tests/test_claims.py), not imported from the engine.
Run from the repository root:  python -m pytest -q tests/test_directional.py
"""
import itertools
import os
import random
import sys

import numpy as np
import sympy as sp

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(ROOT, "directional"))
import axis_timelines as at  # noqa: E402
import clock_rate_vs_sigma as cs  # noqa: E402


# ------------------------------------------------------------ axis timelines
def test_closed_walk_counts_on_z3():
    assert [at.closed_count(n) for n in range(9)] == [1, 0, 6, 0, 90, 0, 1860, 0, 44730]


def test_closed_count_agrees_with_brute_force():
    for n in range(0, 7):
        assert at.closed_count(n) == at.brute_flux_counts(n)[0]


def test_flux_squared_has_the_parity_of_n():
    for n in range(1, 8):
        assert all(f2 % 2 == n % 2 for f2 in at.brute_flux_counts(n))


def test_stopping_rule_is_zero_or_two_for_even_n_and_one_or_three_for_odd_n():
    for n in range(1, 8):
        allowed = {f2 for f2 in at.brute_flux_counts(n) if f2 <= 3}
        assert allowed == ({0, 2} if n % 2 == 0 else {1, 3}) or (n == 1 and allowed == {1})


def test_solver_cost_as_written_adds_the_old_step_where_it_should_remove_it():
    rng = random.Random(3)
    differs = 0
    for _ in range(300):
        tl = [rng.choice(at.DIRS) for _ in range(8)]
        t = rng.randrange(8)
        ax = rng.choice(at.DIRS)
        written, true = at.solver_cost_matches_flux_after_edit(tl, t, ax)
        f = at.flux(tl)
        old, new = at.VEC[tl[t]], at.VEC[ax]
        assert sum((-old[i] + new[i] + f[i]) ** 2 for i in range(3)) == true   # the corrected cost is the true |flux|^2
        differs += written != true
    assert differs > 100


def test_sign_corrected_contract_converges_more_often():
    rows = at.convergence_table(n=12, trials=500, seed=1)
    assert rows["corrected"][0] > rows["as written"][0] + 0.1
    assert rows["corrected"][0] > 0.8


# ------------------------------------------------------------ clock rate
def test_qgd_recipe_clock_is_sqrt_one_minus_sigma_squared():
    """Kerr-Schild g = eta + sigma_t^2 k k^T with k = (1, 1, 0, 0) in (v, r, th, ph): the static clock is sqrt(-g_vv)."""
    s = sp.symbols("sigma_t", positive=True)
    k = sp.Matrix([1, 1, 0, 0])
    eta = sp.diag(-1, 1, 1, 1)
    g = eta + s ** 2 * k * k.T
    assert sp.simplify(g[0, 0] + (1 - s ** 2)) == 0
    # sigma_t^2 = r_s / r  ->  numerical clock
    u = np.array([1.5, 2.0, 10.0])
    assert np.allclose(cs.clock_qgd(u), np.sqrt(1 - 1 / u))


def test_density_that_reproduces_the_qgd_clock_exists_and_diverges_at_the_horizon():
    u = np.array([1.001, 1.5, 3.0, 100.0])
    assert np.allclose(cs.clock_toy(cs.density_required(u)), cs.clock_qgd(u))
    assert cs.density_required(1.0000001) > 1e3


def test_toy_and_qgd_clocks_agree_to_first_order_and_differ_at_second():
    gr, toy, diff = cs.series_difference()
    assert sp.expand(gr - toy) != 0
    assert sp.simplify(sp.diff(diff, cs.x, 0).subs(cs.x, 0)) == 0
    assert sp.simplify(sp.diff(diff, cs.x, 1).subs(cs.x, 0)) == 0
    assert sp.simplify(sp.diff(diff, cs.x, 2).subs(cs.x, 0) - 3) == 0          # difference = 3 x^2 / 2 + ...


def test_lorentzian_well_tail_is_inverse_square_and_the_required_tail_is_inverse_first_power():
    u = sp.symbols("u", positive=True)
    assert cs.tail_exponent(cs.w / (1 + cs.d ** 2), cs.d) == -2
    req = sp.series((1 - 1 / u) ** sp.Rational(-1, 2) - 1, u, sp.oo, 3).removeO()
    assert cs.tail_exponent(req, u) == -1


def test_toy_density_is_bounded_so_its_clock_has_a_floor():
    w = 4.0
    d = np.arange(0, 1000)
    rho = 1.0 + w / (1.0 + d ** 2)
    assert rho.max() == 1.0 + w
    assert cs.clock_toy(rho).min() >= 1.0 / (1.0 + w) - 1e-12 > 0
    assert cs.clock_qgd(1.0) == 0.0


# ------------------------------------------------------------ gravitational relaxation
def test_uniform_density_relaxation_erases_the_forced_event():
    rng = random.Random(7)
    for _ in range(50):
        tl = [rng.choice(at.DIRS) for _ in range(10)]
        tl[0] = "+X"
        out = at.gravity_relax(tl, [1.0] * 10, 2, "+Z", rng)
        assert out[2] != "+Z"
        assert all(a == out[0] for a in out[1:9])            # steps 1..n-2 copy step 0
