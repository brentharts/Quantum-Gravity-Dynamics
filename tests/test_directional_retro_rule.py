"""Tests for directional/retro_rule.py: a one-parameter consistency rule that pays the price of retro-dependence in the hexagonal toy."""
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "directional"))

import retro_rule as R  # noqa: E402

SIX = list(R.V)
TWELVE = list(R.TWELVE)


def test_kappa_zero_is_the_local_model():
    assert abs(abs(R.best_chsh(R.rule("AB", 0.0), SIX)[0]) - 2.0) < 1e-12
    assert R.delta_used(R.rule("AB", 0.0), SIX) < 1e-12


def test_positive_kappa_violates_the_local_bound_on_the_six():
    for kind in ("AB", "align"):
        for k in (0.1, 0.5, 1.0):
            assert abs(R.best_chsh(R.rule(kind, k), SIX)[0]) > 2.0 + 1e-6


def test_chsh_grows_with_kappa_and_passes_the_quantum_value():
    vals = [abs(R.best_chsh(R.rule("AB", k), SIX)[0]) for k in (0.1, 0.2, 0.5, 1.0, 2.0)]
    assert all(b > a for a, b in zip(vals, vals[1:]))
    assert vals[2] > 2 * math.sqrt(2) and vals[-1] > 3.9
    assert abs(R.kappa_for(2 * math.sqrt(2)) - 0.362) < 0.002


def test_no_signalling_on_the_six_for_both_families():
    for kind in ("AB", "align"):
        for k in (0.3, 1.0, 2.0):
            assert R.signalling(R.rule(kind, k), SIX) < 1e-12


def test_no_signalling_is_a_symmetry_property_a_control_signals():
    assert R.signalling(R.rule("oddcontrol", 0.5), SIX) > 0.1


def test_the_rule_signals_when_zero_outcomes_are_allowed():
    for kind in ("AB", "align"):
        assert R.signalling(R.rule(kind, 0.5), TWELVE) > 0.05


def test_consistent_with_the_lp_price():
    for k in (0.1, 0.2, 0.5, 1.0, 2.0):
        P = R.rule("AB", k)
        S = abs(R.best_chsh(P, SIX)[0])
        d = R.delta_used(P, SIX)
        assert S <= min(4.0, 2 + 12 * d) + 1e-9
