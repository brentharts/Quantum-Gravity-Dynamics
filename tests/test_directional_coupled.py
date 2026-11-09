"""Checks of directional/coupled_network.py: the entangled pair of timelines and the CHSH question.
Run from the repository root:  python -m pytest -q tests/test_directional_coupled.py
"""
import math
import os
import sys

import numpy as np

ROOT = os.path.join(os.path.dirname(__file__), "..")
sys.path.insert(0, os.path.join(ROOT, "directional"))
import coupled_network as cn  # noqa: E402
from axis_timelines import DIRS  # noqa: E402


def test_node_contract_returns_the_reference_for_every_pair():
    assert all(cn.node_contract(c, r) == r for c in DIRS for r in DIRS)


def test_measurement_overwrites_the_past_of_a_and_the_whole_of_b():
    a, b = cn.measure_a(*cn.make_pair(8, "X"), 5, "+Z")
    assert a[:6] == ["+Z"] * 6
    assert b == ["-Z"] * 8


def test_b_history_reveals_the_axis_that_a_measured():
    for f in "XYZ":
        for s in "+-":
            _, b = cn.measure_a(*cn.make_pair(8, "X"), 5, s + f)
            assert cn.decode_setting_from_b(b) == f


def test_six_direction_chsh_is_two_for_the_toy_and_for_the_singlet():
    assert abs(cn.six_direction_max(cn.toy_correlation) - 2.0) < 1e-12
    assert abs(cn.six_direction_max(cn.singlet_correlation) - 2.0) < 1e-12


def test_toy_correlation_equals_the_singlet_on_the_six_directions():
    vs = [np.array(v, dtype=float) for v in cn.VEC.values()]
    for a in vs:
        for b in vs:
            assert abs(cn.toy_correlation(a, b) - cn.singlet_correlation(a, b)) < 1e-12


def test_extended_readout_gives_chsh_four_while_the_singlet_gives_two_root_two():
    ex, ey = np.array([1.0, 0, 0]), np.array([0, 1.0, 0])
    b1, b2 = np.array([1, 1, 0]) / math.sqrt(2), np.array([-1, 1, 0]) / math.sqrt(2)
    assert abs(abs(cn.chsh(cn.toy_correlation, ex, ey, b1, b2)) - 4.0) < 1e-12
    assert abs(abs(cn.chsh(cn.singlet_correlation, ex, ey, b1, b2)) - 2 * math.sqrt(2)) < 1e-12


def test_b_outcomes_are_unbiased_whatever_a_measures():
    """No signalling at the level of outcomes: B's readout averaged over A's fair outcome is zero for every setting pair."""
    ax = {v: k for k, v in cn.VEC.items()}
    rng = np.random.default_rng(0)
    for alpha in cn.VEC.values():
        for _ in range(20):
            beta = rng.normal(size=3)
            beta /= np.linalg.norm(beta)
            outs = []
            for s in (+1, -1):
                _, b = cn.measure_a(*cn.make_pair(8, "X"), 5, ax[tuple(int(s * c) for c in alpha)])
                outs.append(np.sign(np.dot(beta, np.array(cn.VEC[b[-1]], dtype=float))))
            assert sum(outs) == 0
