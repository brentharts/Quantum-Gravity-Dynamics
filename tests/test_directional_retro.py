"""Checks of directional/retro_boundary_model.py (settings as boundary conditions).
Run from the repository root:  python -m pytest -q tests/test_directional_retro.py
"""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "directional"))
import contract_options as co  # noqa: E402
import retro_boundary_model as rb  # noqa: E402


def smax_closed_form(delta):
    return min(4.0, 4.0 / 3.0 + 8.0 * delta)


@pytest.mark.parametrize("delta", [0.0, 0.05, 0.1, 0.3, 0.5])
def test_no_signalling_optimum_follows_the_closed_form(delta):
    got = rb.scan([delta])[delta]
    assert got == pytest.approx(smax_closed_form(delta), abs=1e-6)


def test_zero_dependence_is_the_local_model_of_contract_options():
    assert rb.scan([0.0])[0.0] == pytest.approx(co.chsh_max(co.correlation_R), abs=1e-6)
    assert rb.scan([0.0])[0.0] <= 2.0


def test_threshold_for_a_bell_violation_is_one_twelfth():
    quads = list(rb.representative_quadruples())
    assert max(rb.lp_max_chsh(*q, 1 / 12 - 0.005) for q in quads) < 2.0
    assert max(rb.lp_max_chsh(*q, 1 / 12 + 0.005) for q in quads) > 2.0


def test_weight_rule_signals_for_every_positive_kappa_and_not_at_zero():
    assert rb.signalling_of_rule(0.0) < 1e-12
    assert all(rb.signalling_of_rule(k) > 0.05 for k in (0.25, 1.0, 3.0))


def test_weight_rule_exceeds_two_only_with_large_signalling():
    assert rb.rule_chsh(1.0) < 2.0
    assert rb.rule_chsh(3.0) > 2.0 and rb.signalling_of_rule(3.0) > 0.5
