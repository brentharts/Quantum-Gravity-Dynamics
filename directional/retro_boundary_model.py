#!/usr/bin/env python3
"""
retro_boundary_model.py -- settings as boundary conditions: how much must the settings change the history ensemble?

The earlier policies (contract_options.py) fix the history first and let the settings act afterwards.  Here the history ensemble may depend on BOTH
settings, as it would if the settings were boundary conditions at the two ends of the timelines (a retrocausal reading).  The hidden state is the toy's own:
the origin axis v, one of the six directions, with A holding +v and B holding -v (joint flux zero).  Outcomes are read from it:
    out_A(a, v) = sgn(a . v),     out_B(b, v) = sgn(b . (-v)),
with a in the six directions and b in the 26 non-zero vectors of {-1,0,1}^3.  For each pair of settings (a, b) the model supplies a distribution P_ab(v).

Constraints and measure
  no signalling:  A's outcome distribution given a is the same for every b, and B's given b is the same for every a (checked across the four setting pairs that enter one CHSH value).
  dependence budget delta:  |P_ab(v) - 1/6| <= delta for every v, a, b.  delta = 0 means the ensemble ignores the settings (a local hidden-variable model); delta = 5/6 allows anything.

Results (Exact: each is the optimum of a linear program, computed with scipy's HiGHS solver; the tests check the stated bounds)
  1.  delta = 0 gives the local value (4/3 here, as in contract_options.py, and never above 2).
  2.  The largest CHSH value reachable under no signalling is  S_max(delta) = min(4, 4/3 + 8 delta)  (matched by the LP at nine values of delta, from 0 to 5/6).  It exceeds 2 only when
      delta > 1/12, reaches 2 sqrt 2 at delta = 0.187 and the algebraic maximum 4 at delta = 1/3.  So for this state space and these readouts the settings must move each hidden-state probability
      by more than 1/12 (against the uniform 1/6, a relative change of 50%) before any Bell violation is possible: that is the *price of retro-dependence*.
  3.  The weight rule  P_ab(v) proportional to exp(-kappa * out_A * out_B)  (settings favour the origin axes that make the outcomes anti-correlated) signals for every kappa > 0.
Labels: results 1 and 3 Exact; 2 Numeric (an LP optimum to solver tolerance, over the setting quadruples a = +X, a' in {+Y, -X}, b and b' in the 26).
What this is not: it does not say the toy is quantum, and the ternary outcomes (0 for orthogonal) mean the quantum value 2 sqrt 2 is not the natural reference; 2 and 4 are.

    python3 directional/retro_boundary_model.py
"""
import itertools

import numpy as np
from scipy.optimize import linprog

from axis_timelines import VEC

SIX = [np.array(v, dtype=int) for v in VEC.values()]
TWENTY_SIX = [np.array(v, dtype=int) for v in itertools.product((-1, 0, 1), repeat=3) if any(v)]


def out_a(a, v):
    return int(np.sign(np.dot(a, v)))


def out_b(b, v):
    return int(np.sign(np.dot(b, -v)))


def corr_vector(a, b):
    """c[v] = out_A(a, v) out_B(b, v): E(a, b) = sum_v P_ab(v) c[v]."""
    return np.array([out_a(a, v) * out_b(b, v) for v in SIX], dtype=float)


def lp_max_chsh(a, a2, b, b2, delta):
    """max |E(a,b) - E(a,b') + E(a',b) + E(a',b')| under no signalling and |P - 1/6| <= delta; returns the optimum (float)."""
    pairs = [(a, b), (a, b2), (a2, b), (a2, b2)]
    signs = [+1, -1, +1, +1]
    n = 6
    best = 0.0
    for overall in (+1, -1):
        c = np.zeros(4 * n)
        for k, ((x, y), sg) in enumerate(zip(pairs, signs)):
            c[k * n:(k + 1) * n] = -overall * sg * corr_vector(x, y)
        A_eq, b_eq = [], []
        for k in range(4):                                  # normalisation
            row = np.zeros(4 * n)
            row[k * n:(k + 1) * n] = 1
            A_eq.append(row)
            b_eq.append(1.0)

        def add_equal(k1, k2, f):
            for s in (-1, 0, 1):
                row = np.zeros(4 * n)
                for i, v in enumerate(SIX):
                    if f(v) == s:
                        row[k1 * n + i] += 1
                        row[k2 * n + i] -= 1
                A_eq.append(row)
                b_eq.append(0.0)
        # pairs index: 0:(a,b) 1:(a,b') 2:(a',b) 3:(a',b')
        add_equal(0, 1, lambda v: out_a(a, v))              # A at a: b against b'
        add_equal(2, 3, lambda v: out_a(a2, v))             # A at a'
        add_equal(0, 2, lambda v: out_b(b, v))              # B at b: a against a'
        add_equal(1, 3, lambda v: out_b(b2, v))             # B at b'
        bounds = [(max(0.0, 1 / 6 - delta), min(1.0, 1 / 6 + delta))] * (4 * n)
        res = linprog(c, A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds, method="highs")
        if res.status == 0:
            best = max(best, -res.fun)
    return best


def representative_quadruples():
    ex, ey = SIX[0], SIX[2]
    for a2 in (ey, -ex):
        for b, b2 in itertools.product(TWENTY_SIX, repeat=2):
            yield ex, a2, b, b2


def scan(deltas, quadruples=None):
    """{delta: largest CHSH over the quadruples}."""
    quads = list(quadruples if quadruples is not None else representative_quadruples())
    return {d: max(lp_max_chsh(a, a2, b, b2, d) for a, a2, b, b2 in quads) for d in deltas}


# ---------------------------------------------------------------- a natural weight rule
def weight_rule(a, b, kappa):
    w = np.exp(-kappa * corr_vector(a, b))
    return w / w.sum()


def signalling_of_rule(kappa):
    """Largest change of A's outcome distribution when only b changes (and of B's when only a changes), over the settings."""
    worst = 0.0
    for a in SIX:
        dists = []
        for b in TWENTY_SIX:
            P = weight_rule(a, b, kappa)
            dists.append([sum(P[i] for i, v in enumerate(SIX) if out_a(a, v) == s) for s in (-1, 0, 1)])
        d = np.array(dists)
        worst = max(worst, np.abs(d - d[0]).max())
    return worst


def rule_chsh(kappa):
    """Largest CHSH value of the weight rule over the representative setting quadruples (it signals for kappa > 0)."""
    best = 0.0
    for a, a2, b, b2 in representative_quadruples():
        E = lambda x, y: float(np.dot(weight_rule(x, y, kappa), corr_vector(x, y)))
        best = max(best, abs(E(a, b) - E(a, b2) + E(a2, b) + E(a2, b2)))
    return best


def main():
    import time
    t0 = time.time()
    deltas = [0.0, 0.02, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 5 / 6]
    print("delta  max CHSH under no signalling")
    for d, s in scan(deltas).items():
        print("%.3f   %.4f" % (d, s))
    print("(%.0f s)" % (time.time() - t0))
    for kappa in (0.0, 0.25, 1.0):
        print("weight rule kappa = %.2f: signalling %.4f, max CHSH %.4f" % (kappa, signalling_of_rule(kappa), rule_chsh(kappa)))
    print("closed form  S_max(delta) = min(4, 4/3 + 8 delta): exceeds 2 at delta = 1/12 = %.4f, reaches 2 sqrt2 at %.4f, reaches 4 at 1/3" %
          (1 / 12, (2 * 2 ** 0.5 - 4 / 3) / 8))


if __name__ == "__main__":
    main()
