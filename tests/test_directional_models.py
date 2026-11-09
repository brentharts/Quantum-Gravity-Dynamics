"""Checks of directional/contract_options.py, two_source_clocks.py and update_rate_clock.py.
Run from the repository root:  python -m pytest -q tests/test_directional_models.py
"""
import math
import os
import sys

import numpy as np
import pytest
import sympy as sp

HERE = os.path.dirname(__file__)
sys.path.insert(0, os.path.join(HERE, "..", "directional"))
sys.path.insert(0, HERE)
import contract_options as co  # noqa: E402
import two_source_clocks as ts  # noqa: E402
import update_rate_clock as uc  # noqa: E402


# ---------------------------------------------------------------- contract options
def test_global_conservation_signals_through_b_flux():
    ex, ey = co.SIX[0], co.SIX[2]
    assert co.total_variation(co.b_flux_distribution(ex, 6), co.b_flux_distribution(ey, 6)) > 0.1
    assert co.total_variation(co.b_flux_distribution(ex, 6), co.b_flux_distribution(-ex, 6)) < 1e-12


def test_local_conservation_has_no_correlation_and_no_signalling():
    assert co.chsh_max(co.correlation_L) == 0.0


def test_revelation_obeys_the_local_bound():
    assert co.chsh_max(co.correlation_R) <= 2.0 + 1e-12


def test_copy_policy_reaches_four():
    assert abs(co.chsh_max(co.correlation_C) - 4.0) < 1e-12


def test_global_policy_flux_readout_stays_below_two():
    assert co.chsh_max(co.correlation_G) < 2.0


# ---------------------------------------------------------------- two sources
def test_sigma_sum_clock_differs_from_the_additive_one_at_first_order_in_G():
    ser = ts.first_order_cross_term()
    eps = sp.symbols("epsilon", positive=True)
    coeff = sp.expand(ser).coeff(eps, 2)
    assert sp.simplify(coeff + ts.x1 + ts.x2 + 2 * sp.sqrt(ts.x1) * sp.sqrt(ts.x2)) == 0


def test_additive_density_and_linear_gr_agree_to_first_order():
    x = np.array([1e-4, 3e-4])
    assert abs(ts.clock_toy(list(x)) - ts.clock_linear_gr(list(x))) < 2 * x.sum() ** 2
    assert abs(ts.clock_sigma_sum(list(x)) - ts.clock_linear_gr(list(x))) > 0.5 * 2 * math.sqrt(x[0] * x[1])


def test_far_field_mass_of_the_sigma_sum_is_the_square_of_the_sum_of_square_roots():
    ratio = ts.far_field_mass_ratio()
    assert sp.simplify(ratio - (sp.sqrt(ts.M1) + sp.sqrt(ts.M2)) ** 2 / (ts.M1 + ts.M2)) == 0
    assert ratio.subs({ts.M1: 1, ts.M2: 1}) == 2


def test_cross_term_is_not_harmonic_but_a_single_source_is():
    lap_cross, lap_single = ts.laplacian_cross_term_at()
    assert abs(lap_single) < 1e-12
    assert abs(lap_cross) > 1e-2


def test_sigma_sum_ricci_residual_is_first_order_in_G_not_second():
    """Same helpers and points as tests/test_gaps.py; halving both masses halves the sigma-sum residual and quarters the linear Kerr-Schild one."""
    import test_gaps as tg

    def resid(m, sigma_sum):
        parts = [tg.ks_one_form_term(m, (-3.0, 0, 0)), tg.ks_one_form_term(m, (3.0, 0, 0))]
        if sigma_sum:
            s = sum((sp.sqrt(h) * k for h, k in parts), sp.zeros(4, 1))
            g = tg.ETA + s * s.T
        else:
            g = tg.ETA + sum((h * k * k.T for h, k in parts), sp.zeros(4, 4))
        R, _ = tg.ricci_at(g, tg.P)
        return np.max(np.abs(R))

    assert resid(0.005, True) / resid(0.01, True) == pytest.approx(0.5, rel=0.05)
    assert resid(0.005, False) / resid(0.01, False) == pytest.approx(0.25, rel=0.05)


# ---------------------------------------------------------------- update-rate clock
def test_tick_count_matches_its_binomial_expectation():
    rng = np.random.default_rng(11)
    for rho in (1.0, 1.2247, 3.0):
        n = uc.simulate_ticks(rho, 400_000, rng)
        mean, sd = uc.tick_statistics(rho, 400_000)
        assert abs(n - mean) <= 5 * max(sd, 1.0)


def test_received_over_emitted_ticks_is_the_density_ratio():
    rng = np.random.default_rng(5)
    rho1, rho2 = uc.rho_qgd(3.0), uc.rho_qgd(10.0)
    got = np.mean([uc.received_ticks(rho1, rho2, 200, rng) for _ in range(2000)]) / 200
    assert got == pytest.approx(rho1 / rho2, rel=0.01)


def test_qgd_matching_density_gives_the_schwarzschild_redshift():
    assert uc.rho_qgd(3.0) / uc.rho_qgd(10.0) == pytest.approx(uc.redshift_prediction(3.0, 10.0), rel=1e-12)
