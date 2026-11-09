"""Tests for directional/hex_directions.py: the six-direction toy with the hexagonal (red A2) direction set."""
import math
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "directional"))

import hex_directions as H  # noqa: E402


def test_closed_timeline_counts_hex_and_cubic():
    assert [H.closed_count_hex(n) for n in range(9)] == [1, 0, 6, 12, 90, 360, 2040, 10080, 54810]
    assert [H.closed_count_cubic(n) for n in range(9)] == [1, 0, 6, 0, 90, 0, 1860, 0, 44730]


def test_hex_count_agrees_with_brute_force():
    for n in range(0, 7):
        assert H.brute_closed_hex(n) == H.closed_count_hex(n)


def test_hex_has_closed_three_step_timelines_and_cubic_has_none():
    assert H.minimal_closed_words(3, H.HINT) == 12
    assert H.minimal_closed_words(3, H.CUBIC_INT) == 0
    assert H.hflux_is_zero(["+A", "+C", "-B"])          # e_0 + e_2 - e_1 = 0


def test_copy_contract_signals_the_axis_family_in_the_hexagonal_set():
    for d in H.HDIRS:
        a, b = H.hex_measure_a(["+A"] + ["+B"] * 7, ["-A"] + ["+B"] * 7, 5, d)
        assert b[-1] == H.hinvert(d)
        assert H.hex_setting_family_from_b(b) == d[1]
    families = {H.hex_setting_family_from_b(H.hex_measure_a(["+A"] + ["+B"] * 7, ["-A"] + ["+B"] * 7, 5, d)[1]) for d in H.HDIRS}
    assert len(families) == 3                              # log2 3 bits


def test_no_pair_of_hexagonal_directions_is_orthogonal():
    assert all(H.hdot(a, b) != 0 for a in H.HDIRS for b in H.HDIRS)


def test_toy_correlation_is_pm_one_on_the_six_and_reaches_chsh_four_while_the_singlet_does_not():
    assert all(abs(H.toy_E(a, b)) == 1 for a in H.HDIRS for b in H.HDIRS)
    assert H.chsh_max_six(H.toy_E) == 4.0
    assert H.chsh_max_six(H.singlet_E) == pytest.approx(2.5)
    assert H.chsh_max_six(H.singlet_E) < 2 * math.sqrt(2)


@pytest.mark.parametrize("delta", [0.0, 0.02, 0.05, 0.1, 0.15, 0.1667, 0.2])
def test_extremal_quadruple_attains_two_plus_twelve_delta(delta):
    a, a2, b, b2 = H.SIX_V[0], H.SIX_V[2], H.SIX_V[3], H.SIX_V[5]        # e_0, e_1 (60 deg), -e_0 (180), -e_2 (300)
    assert H.lp_max_chsh_h(a, a2, b, b2, delta) == pytest.approx(min(4.0, 2 + 12 * delta), abs=1e-6)


@pytest.mark.parametrize("delta", [0.0, 0.05, 0.1])
def test_no_quadruple_of_the_scan_exceeds_the_formula(delta):
    got = H.scan_hex([delta])[delta]            # a = e_0, a' in {e_1, e_2, -e_0}, b and b' among the 12 directions
    assert got == pytest.approx(min(4.0, 2 + 12 * delta), abs=1e-6)


def test_any_dependence_gives_a_violation_in_the_hexagonal_set():
    a, a2, b, b2 = H.SIX_V[0], H.SIX_V[2], H.SIX_V[3], H.SIX_V[5]
    assert H.lp_max_chsh_h(a, a2, b, b2, 0.01) > 2.0
