#!/usr/bin/env python3
"""
coupled_network.py -- the entangled pair of six-direction timelines, and whether a Bell (CHSH) test can say anything about it.

Port of CoupledParticleNetwork (spectre-monotile-py/research/toymodel_v1) without printing.  Two timelines A and B start anti-correlated
(A[0] = +axis, B[0] = -axis).  A measurement of A at step m forces A[m] to the observed axis, copies it backward to A[0], sets B[0] to the
inverse of A[0], and copies it forward along B.

Results (Exact unless stated; every one is a test in tests/test_directional.py)
-------
Exact 1.  The node contract `_resolve_node_contract(current, reference)` returns `reference` for every pair (36 pairs).  The contract does
          not use the current state; "contracting" a node is copying its neighbour.  So after A is measured along an axis, all of A[0..m]
          equal the observed axis and all of B equal its inverse.
Exact 2.  B's timeline after the measurement determines the axis family (X, Y or Z) that A measured, with certainty.  An observer of B's history
          would learn A's *setting* (log2 3 bits) instantly: as written, the model signals.
Exact 3.  If both sides read out only along the six directions, the toy's correlation is E(a, b) = -a.b, the same as the quantum singlet's,
          and the largest CHSH value over all settings from the six directions is 2 for both.  So no CHSH experiment restricted to the six
          directions can tell the toy from quantum mechanics or from a local model.
Exact 4.  The toy has no rule for what a measurement of B along another direction returns.  Extending it by the simplest rule (the outcome is the
          sign of the projection of B's vector on the chosen direction; 0 for orthogonal) gives E(a, b) = -sgn(a.b) for A restricted to the six
          directions.  With a = +X, a' = +Y, b = (+X+Y)/sqrt2, b' = (-X+Y)/sqrt2, |S| = 4, above both the local bound 2 and the quantum bound
          2 sqrt2 (the quantum singlet gives 2 sqrt2 for the same settings).  Each party's own outcomes stay uniform, so the *outcome* statistics do not
          signal; Exact 2 is the signalling, and it is in the history.
Hypothesis.  That a different contract (one that uses the current state, or reveals B's history only through a Born-like rule) removes the
          signalling and reaches 2 sqrt2: not built here.

    python3 directional/coupled_network.py
"""
import itertools
import math

import numpy as np

from axis_timelines import DIRS, VEC


def invert(axis):
    return ("-" if axis[0] == "+" else "+") + axis[1]


def node_contract(current, reference):
    """As written: the direction closest to the reference vector (current is not used)."""
    rv = VEC[reference]
    return min(DIRS, key=lambda ax: sum((VEC[ax][i] - rv[i]) ** 2 for i in range(3)))


def make_pair(length=8, axis="X", fill="+Y"):
    return [("+" + axis)] + [fill] * (length - 1), [("-" + axis)] + [fill] * (length - 1)


def measure_a(tl_a, tl_b, m, observed):
    """CoupledParticleNetwork.measure_particle_a: returns the new (A, B)."""
    a, b = list(tl_a), list(tl_b)
    a[m] = observed
    for t in range(m - 1, -1, -1):
        a[t] = node_contract(a[t], a[t + 1])
    b[0] = invert(a[0])
    for t in range(1, len(b)):
        b[t] = node_contract(b[t], b[t - 1])
    return a, b


def decode_setting_from_b(tl_b):
    """What an observer of B's history reads: the axis family of its last step."""
    return tl_b[-1][1]


def toy_correlation(alpha, beta):
    """E(alpha, beta) for A measured along the unit vector alpha (one of the six) and B read along the vector beta,
    with outcome(A) = s uniform in {+1, -1} and outcome(B) = sgn(beta . B_last)."""
    ax = {v: k for k, v in VEC.items()}
    total = 0.0
    for s in (+1, -1):
        observed = ax[tuple(int(s * c) for c in alpha)]
        _, b = measure_a(*make_pair(8, "X"), 5, observed)
        bv = np.array(VEC[b[-1]], dtype=float)
        out_b = np.sign(np.dot(beta, bv))
        total += s * out_b
    return total / 2.0


def singlet_correlation(alpha, beta):
    return -float(np.dot(alpha, beta))


def chsh(E, a, a2, b, b2):
    return E(a, b) - E(a, b2) + E(a2, b) + E(a2, b2)


def six_direction_max(E):
    vs = [np.array(v, dtype=float) for v in VEC.values()]
    return max(abs(chsh(E, a, a2, b, b2)) for a, a2, b, b2 in itertools.product(vs, repeat=4))


def main():
    print("Exact 1: node_contract(current, reference) == reference for all 36 pairs:",
          all(node_contract(c, r) == r for c in DIRS for r in DIRS))
    a, b = measure_a(*make_pair(8, "X"), 5, "+Z")
    print("after A measured +Z at step 5:\n  A =", a, "\n  B =", b)
    print("Exact 2: B decodes the axis family A measured:", all(decode_setting_from_b(measure_a(*make_pair(8, "X"), 5, s + f)[1]) == f
                                                          for f in "XYZ" for s in "+-"))
    print("Exact 3: max CHSH over the six directions, toy %.3f, singlet %.3f" % (six_direction_max(toy_correlation), six_direction_max(singlet_correlation)))
    ex, ey = np.array([1.0, 0, 0]), np.array([0, 1.0, 0])
    b1, b2 = np.array([1, 1, 0]) / math.sqrt(2), np.array([-1, 1, 0]) / math.sqrt(2)
    print("Exact 4: a = X, a' = Y, b = (X+Y)/sqrt2, b' = (-X+Y)/sqrt2:  toy S = %.3f, singlet S = %.3f (2 sqrt2 = %.3f)"
          % (chsh(toy_correlation, ex, ey, b1, b2), chsh(singlet_correlation, ex, ey, b1, b2), 2 * math.sqrt(2)))


if __name__ == "__main__":
    main()
