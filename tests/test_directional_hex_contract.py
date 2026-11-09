"""Tests for directional/hex_contract.py: shock repair on the hexagonal timeline, with and without a 3-cycle reward."""
import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "directional"))

import hex_contract as C  # noqa: E402


def test_cubic_even_length_always_repairs_with_one_other_edit():
    assert dict(C.repair_distribution(C.CUB, 6)) == {1: 1860 * 6 * 5}


def test_cubic_odd_length_has_no_closed_timeline_so_nothing_to_repair():
    assert len(C.closed_words(C.CUB, 5)) == 0
    assert dict(C.repair_distribution(C.CUB, 5)) == {}


def test_hexagonal_repair_needs_one_or_two_other_edits():
    d5 = C.repair_distribution(C.HEX, 5)
    d6 = C.repair_distribution(C.HEX, 6)
    assert dict(d5) == {1: 7560, 2: 1440}
    assert dict(d6) == {1: 57960, 2: 3240}
    assert sum(d5.values()) == 360 * 5 * 5 and sum(d6.values()) == 2040 * 6 * 5


def test_hexagonal_closed_timelines_exist_for_odd_lengths():
    assert [len(C.closed_words(C.HEX, n)) for n in (3, 5)] == [12, 360]


@pytest.mark.parametrize("kind,mu,closes,edits,tri", [
    ("cub", 0, 1.0, 1.133, 0.0), ("hex", 0, 0.922, 1.215, 0.225), ("hex", 0.5, 0.920, 1.494, 0.814), ("hex", 1.5, 0.731, 2.411, 2.158)])
def test_greedy_solver_rows_n6(kind, mu, closes, edits, tri):
    r = C.run(kind, 6, mu)
    assert r["closes"] == pytest.approx(closes, abs=0.002)
    assert r["edits"] == pytest.approx(edits, abs=0.002)
    assert r["T"] == pytest.approx(tri, abs=0.002)


def test_reward_raises_triples_and_edits_and_strong_reward_costs_closure():
    r0, r1, r2 = C.run("hex", 8, 0), C.run("hex", 8, 0.5), C.run("hex", 8, 1.5)
    assert r0["T"] < r1["T"] < r2["T"]
    assert r0["edits"] < r1["edits"] < r2["edits"]
    assert r0["reach"] < r2["reach"]
    assert r2["closes"] < 0.6 < r0["closes"]
