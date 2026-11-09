"""Tests for directional/hex_e8.py: the six hexagonal directions as roots of E8."""
import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "directional"))

import hex_directions as H  # noqa: E402
import hex_e8 as E  # noqa: E402


def test_e8_has_240_roots_of_norm_two():
    R = E.roots()
    assert R.shape == (240, 8) and np.allclose((R ** 2).sum(1), 2)


def test_the_six_hexagonal_vectors_are_roots_and_the_map_is_an_isometry_up_to_two():
    A = E.hex_roots()
    keys = {tuple(np.round(r, 6) + 0.0) for r in E.roots()}
    assert all(tuple(np.round(a, 6) + 0.0) in keys for a in A)
    G = A @ A.T
    names = list(E.HEX_PQ)
    want = np.array([[2 * H.hdot(a, b) for b in names] for a in names])
    assert np.allclose(G, want)


def test_orthogonal_complement_has_72_roots():
    assert len(E.orthogonal_roots()) == 72


def test_sign_patterns_of_the_240_roots_on_the_six():
    pats = E.distinct_readouts()
    assert len(pats) == 13
    assert sum(0 in p for p in pats) == 7


@pytest.mark.parametrize("delta", [0.05, 0.2])
def test_lp_with_e8_readouts_follows_two_plus_twelve_delta(delta):
    assert E.scan_e8([delta])[delta] == pytest.approx(min(4.0, 2 + 12 * delta), abs=1e-6)
