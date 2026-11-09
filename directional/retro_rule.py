#!/usr/bin/env python3
"""
retro_rule.py -- a consistency rule that pays the price of retro-dependence in the hexagonal toy, and what it does and does not show.

Why.  retro_boundary_model.py / hex_directions.py give the PRICE of a Bell violation (the largest CHSH value under no signalling is
min(4, 2 + 12 delta), delta = how far each hidden-state probability may move from 1/6) but no rule that pays it: the LEDGER row "A consistency rule
(not a free choice of P_ab)" was UNTESTED.  Here a rule is a map (a, b) -> P(v | a, b), v = origin direction among the six.  Two families, each with one
parameter kappa >= 0 (kappa = 0 is the uniform, local model):

    AB     P(v | a, b) proportional to exp(kappa * A(a, v) * B(b, v))     A, B = the +-1/0 sign read-outs: a Boltzmann weight on the
                                                                          agreement of the two ends of the history, normalised per (a, b)
    align  P(v | a, b) proportional to exp(kappa * (a.v) (b.v))           the same with the inner products instead of the signs

This is the "boundary terms at both ends and a normalisation over the whole history" candidate of docs/DIRECTIONAL_AXES.md: the dependence on the two settings
is not put in as a table of probabilities but as one coupling between the two ends.

Results (tests in tests/test_directional_retro_rule.py)
Exact R1.  kappa = 0 is the local model (best S = 2 exactly).  For kappa > 0 and read-outs among the six directions, S exceeds 2.
Exact R2.  No signalling holds exactly for read-outs among the six hexagonal directions.  This is a consequence of the symmetry v -> -v of the weights, not of
           the particular rule: with no zero outcomes each party's marginal is exactly uniform on +-1.  A control rule that breaks the symmetry
           (weight exp(kappa (a.v + b.v))) signals.
Exact R3.  With read-outs among the twelve directions (where a . v = 0 can occur) the same rules SIGNAL (the probability of the zero outcome depends on the
           other setting).  So the rule pays the price only with zero-free read-outs, which restricts the readout set the same way the cubic zeros did.
Numeric R4. CHSH against kappa (AB, best over the settings among the six): see the table printed by main(); S exceeds the quantum value 2 sqrt2 at kappa = 0.362
           and approaches 4 for large kappa (3.62 at kappa = 1, 3.95 at 2).  A value above 2 sqrt2 is not a quantum value: the rule is a no-signalling box, not a model of the singlet.
Exact R5.  Every value is consistent with the LP bound S <= min(4, 2 + 12 delta(kappa)), delta(kappa) = max |P - 1/6| over the settings used.
What this is not.  The weight is a CHOICE (one parameter kappa and a family), not derived from the path action of the formal section; it puts the dependence
of the hidden state on both settings in by fiat, with a per-pair normalisation Z(a, b); and 'no signalling' here is a symmetry property.  It shows that the
price can be paid by a one-parameter rule, not that physics pays it.
"""
import itertools
import math
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import hex_directions as H  # noqa: E402

V = H.SIX_V
TWELVE = H.TWELVE


def _sg(x):
    return np.where(np.abs(x) < 1e-9, 0.0, np.sign(x))


def weights(kind, kappa, a, b):
    pa, pb = V @ a, V @ b
    if kind == "AB":
        w = np.exp(kappa * _sg(pa) * _sg(pb))
    elif kind == "align":
        w = np.exp(kappa * pa * pb)
    elif kind == "oddcontrol":                      # breaks v -> -v: a control
        w = np.exp(kappa * (pa + pb))
    else:
        raise ValueError(kind)
    return w / w.sum()


def rule(kind, kappa):
    return lambda a, b: weights(kind, kappa, a, b)


def E(P, a, b):
    return float(P(a, b) @ (_sg(V @ a) * _sg(V @ b)))


def chsh(P, a, a2, b, b2):
    return E(P, a, b) + E(P, a, b2) + E(P, a2, b) - E(P, a2, b2)


def best_chsh(P, readouts):
    best, arg = 0.0, None
    for i, j, k, l in itertools.product(range(len(readouts)), repeat=4):
        S = chsh(P, readouts[i], readouts[j], readouts[k], readouts[l])
        if abs(S) > abs(best) + 1e-12:
            best, arg = S, (i, j, k, l)
    return best, arg


def signalling(P, readouts):
    """largest change of one party's full outcome distribution (outcomes -1, 0, +1) when the other party's setting changes."""
    s = 0.0
    for a in readouts:
        ms = np.array([[P(a, b)[_sg(V @ a) == t].sum() for t in (-1, 0, 1)] for b in readouts])
        s = max(s, np.abs(ms - ms[0]).max())
    for b in readouts:
        ms = np.array([[P(a, b)[_sg(V @ b) == t].sum() for t in (-1, 0, 1)] for a in readouts])
        s = max(s, np.abs(ms - ms[0]).max())
    return s


def delta_used(P, readouts):
    return max(np.abs(P(a, b) - 1 / 6).max() for a in readouts for b in readouts)


def kappa_for(target, kind="AB", lo=0.0, hi=10.0):
    """smallest kappa (bisection) with best CHSH on the six equal to target."""
    six = list(V)
    for _ in range(40):
        mid = (lo + hi) / 2
        if abs(best_chsh(rule(kind, mid), six)[0]) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def main():
    six = list(V)
    print("AB rule, read-outs among the six: kappa, delta, best |S|, LP bound min(4, 2 + 12 delta), signalling (six), signalling (twelve)")
    for k in (0.0, 0.1, 0.2, 0.34, 0.5, 1.0, 2.0, 5.0):
        P = rule("AB", k)
        S = abs(best_chsh(P, six)[0])
        d = delta_used(P, six)
        print("  %.2f  %.3f  %.3f  %.3f  %.3f  %.3f" % (k, d, S, min(4, 2 + 12 * d), signalling(P, six), signalling(P, list(TWELVE))))
    print("kappa at which S reaches 2 sqrt2 (AB): %.3f" % kappa_for(2 * math.sqrt(2)))


if __name__ == "__main__":
    main()
