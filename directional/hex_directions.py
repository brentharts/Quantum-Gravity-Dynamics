#!/usr/bin/env python3
"""
hex_directions.py -- the same toy with the six HEXAGONAL directions (the red A2 roots of the Spectre lift) instead of +-X, +-Y, +-Z.

Why.  quantgrav.tex says each direction maps to an E8 root, but the embedding of spectre_in_e8 is of A2 + A2, whose red part is six signed directions
in a plane, 60 degrees apart (e_0, e_1, e_2 and their negatives, e_0 - e_1 + e_2 = 0).  This module reruns the questions of axis_timelines.py,
coupled_network.py and retro_boundary_model.py with that set and records what changes.  Nothing here uses E8 beyond the direction set, and nothing
says the toy is quantum or gravitational.

Results (Exact unless stated; all are tests in tests/test_directional_hex.py)
Exact H1.  Closed timelines (net flux zero) of n steps: 1, 0, 6, 12, 90, 360, 2040, 10080 for n = 0..7 (closed walks on the triangular lattice), against
           1, 0, 6, 0, 90, 0, 1860, 0 for the cubic set.  Brute force for n <= 7 agrees with the generating count.
Exact H2.  Parity.  The cubic set is bipartite: |flux|^2 has the parity of n, so a timeline of odd length can never close.  The hexagonal set has closed
           3-step timelines (the 12 ordered triples (e_0, e_2, -e_1) and relatives, e.g. e_0 + e_2 + (-e_1) = 0), so closure has no parity.  The minimal
           closed words (no proper cyclic run closes; ordered count) are 6 of length 2, 12 of length 3 and 24 of length 4 for the hexagonal set, against 6, 0, 24
           for the cubic set: the 3-cycle is new.
Exact H3.  The copy contract signals here too: after A is measured along any of the six, B's last step reveals which of the 3 axes (log2 3 bits) was measured.
Exact H4.  The correlation read out with the sign rule is E(a, b) = -sgn(a . b).  Among the six hexagonal directions no pair is orthogonal, so E = +-1
           for every pair (it is 0 for orthogonal pairs in the cubic case).  Consequence: with all four settings among the six directions the toy
           reaches |CHSH| = 4, while the singlet value -a . b = -cos(angle) gives at most 5/2 on the same six directions, below 2 sqrt2.  The cubic toy
           and the singlet both give 2 on the six cubic directions, so a Bell test restricted to the six directions could not tell them apart there.  Here it can.
Numeric H5.  Price of retro-dependence (the LP of retro_boundary_model.py): hidden state = origin direction among the six, outcomes sgn(a . v), a among
           the six, b among the six or among the 12 directions 30 degrees apart (the same optimum for both); the largest CHSH under no signalling is
               S_max(delta) = min(4, 2 + 12 delta)         (LP optimum, checked at eight values of delta)
           so the local value is already 2 and ANY dependence delta > 0 gives a Bell violation; 2 sqrt2 is reached at delta = (2 sqrt2 - 2)/12 = 0.069 and 4 at 1/6.
           The cubic result (threshold delta > 1/12, S_max = min(4, 4/3 + 8 delta)) therefore depends on the zero outcomes of the cubic read-out (a . v = 0
           for orthogonal pairs), which lower its local value to 4/3; it is not a property of "six directions" as such.
What this is not: a derivation of any of the directions from E8, a physical model, or a statement about the continuum; and 'sgn' with 12 readout
directions is an assumption, as in the cubic case.

    python3 directional/hex_directions.py
"""
import itertools
import math
from collections import Counter

import numpy as np
from scipy.optimize import linprog

HDIRS = ("+A", "-A", "+B", "-B", "+C", "-C")
_ANG = {"+A": 0, "-A": 180, "+B": 60, "-B": 240, "+C": 120, "-C": 300}
HVEC = {k: (round(math.cos(math.radians(a)), 12), round(math.sin(math.radians(a)), 12)) for k, a in _ANG.items()}
# exact integer coordinates in the basis (e_0, e_1): a vector x e_0 + y e_1 has Gram matrix [[1, 1/2], [1/2, 1]]
HINT = {"+A": (1, 0), "-A": (-1, 0), "+B": (0, 1), "-B": (0, -1), "+C": (-1, 1), "-C": (1, -1)}


def hflux(timeline):
    return tuple(sum(HINT[a][i] for a in timeline) for i in range(2))


def hflux_is_zero(timeline):
    return hflux(timeline) == (0, 0)


# ---------------------------------------------------------------- H1, H2: closed timelines
def closed_count_hex(n):
    """Closed walks of n steps on the triangular lattice (six neighbours), by dynamic programming on the integer coordinates."""
    cur = Counter({(0, 0): 1})
    for _ in range(n):
        nxt = Counter()
        for (x, y), c in cur.items():
            for dx, dy in HINT.values():
                nxt[(x + dx, y + dy)] += c
        cur = nxt
    return cur[(0, 0)]


def closed_count_cubic(n):
    cur = Counter({(0, 0, 0): 1})
    steps = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
    for _ in range(n):
        nxt = Counter()
        for p, c in cur.items():
            for d in steps:
                nxt[(p[0] + d[0], p[1] + d[1], p[2] + d[2])] += c
        cur = nxt
    return cur[(0, 0, 0)]


def brute_closed_hex(n):
    return sum(hflux_is_zero(w) for w in itertools.product(HDIRS, repeat=n))


def minimal_closed_words(n, vec):
    """Number of ordered closed words of length n with no proper cyclic contiguous run closing (vec: dict name -> integer tuple)."""
    names = list(vec)
    cnt = 0
    for w in itertools.product(names, repeat=n):
        S = [vec[a] for a in w]
        if any(sum(s[i] for s in S) for i in range(len(S[0]))):
            continue
        ok = True
        for i in range(n):
            tot = [0] * len(S[0])
            for ln in range(1, n):
                s = S[(i + ln - 1) % n]
                tot = [tot[k] + s[k] for k in range(len(tot))]
                if not any(tot):
                    ok = False
                    break
            if not ok:
                break
        cnt += ok
    return cnt


CUBIC_INT = {"+X": (1, 0, 0), "-X": (-1, 0, 0), "+Y": (0, 1, 0), "-Y": (0, -1, 0), "+Z": (0, 0, 1), "-Z": (0, 0, -1)}


# ---------------------------------------------------------------- H3: the copy contract
def hinvert(a):
    return ("-" if a[0] == "+" else "+") + a[1]


def hex_measure_a(tl_a, tl_b, m, observed):
    """The copy contract of coupled_network.py: a[t] copies a[t+1] back to 0, b[0] = inverse of a[0], b[t] copies b[t-1]."""
    a, b = list(tl_a), list(tl_b)
    a[m] = observed
    for t in range(m - 1, -1, -1):
        a[t] = a[t + 1]
    b[0] = hinvert(a[0])
    for t in range(1, len(b)):
        b[t] = b[t - 1]
    return a, b


def hex_setting_family_from_b(tl_b):
    return tl_b[-1][1]


# ---------------------------------------------------------------- H4: correlations and CHSH on the six directions
def hdot(a, b):
    """Exact inner product of two of the six directions: 1, 1/2, -1/2, -1 (as a float of those values)."""
    x, y = HINT[a], HINT[b]
    return x[0] * y[0] + x[1] * y[1] + 0.5 * (x[0] * y[1] + x[1] * y[0])


def sgn(x):
    return (x > 1e-9) - (x < -1e-9)


def toy_E(a, b):
    """Copy contract + sign read-out: A measured along a with outcome s, B holds -s a; out_B = sgn(b . (-s a)); E = mean of s * out_B = -sgn(a . b)."""
    return -sgn(hdot(a, b))


def singlet_E(a, b):
    return -hdot(a, b)


def chsh_max_six(E):
    best = 0.0
    for a, a2, b, b2 in itertools.product(HDIRS, repeat=4):
        best = max(best, abs(E(a, b) - E(a, b2) + E(a2, b) + E(a2, b2)))
    return best


# ---------------------------------------------------------------- H5: the LP
SIX_V = [np.array(HVEC[k]) for k in HDIRS]
TWELVE = [np.array((math.cos(math.radians(30 * k)), math.sin(math.radians(30 * k)))) for k in range(12)]


def _sg(x):
    return int((x > 1e-9) - (x < -1e-9))


def out_a_h(a, v):
    return _sg(float(np.dot(a, v)))


def out_b_h(b, v):
    return _sg(float(np.dot(b, -v)))


def corr_vector_h(a, b):
    return np.array([out_a_h(a, v) * out_b_h(b, v) for v in SIX_V], dtype=float)


def lp_max_chsh_h(a, a2, b, b2, delta):
    pairs = [(a, b), (a, b2), (a2, b), (a2, b2)]
    signs = [+1, -1, +1, +1]
    n = 6
    best = 0.0
    for overall in (+1, -1):
        c = np.zeros(4 * n)
        for k, ((x, y), sg) in enumerate(zip(pairs, signs)):
            c[k * n:(k + 1) * n] = -overall * sg * corr_vector_h(x, y)
        A_eq, b_eq = [], []
        for k in range(4):
            row = np.zeros(4 * n)
            row[k * n:(k + 1) * n] = 1
            A_eq.append(row)
            b_eq.append(1.0)

        def add_equal(k1, k2, f):
            for s in (-1, 0, 1):
                row = np.zeros(4 * n)
                for i, v in enumerate(SIX_V):
                    if f(v) == s:
                        row[k1 * n + i] += 1
                        row[k2 * n + i] -= 1
                A_eq.append(row)
                b_eq.append(0.0)
        add_equal(0, 1, lambda v: out_a_h(a, v))
        add_equal(2, 3, lambda v: out_a_h(a2, v))
        add_equal(0, 2, lambda v: out_b_h(b, v))
        add_equal(1, 3, lambda v: out_b_h(b2, v))
        bounds = [(max(0.0, 1 / 6 - delta), min(1.0, 1 / 6 + delta))] * (4 * n)
        res = linprog(c, A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds, method="highs")
        if res.status == 0:
            best = max(best, -res.fun)
    return best


def quadruples_h(readouts=None):
    """a = e_0, a' in {e_1 (60 deg), e_2 (120 deg), -e_0}; b, b' among the readout directions (default: the 12)."""
    readouts = TWELVE if readouts is None else readouts
    e0, e1, e2 = SIX_V[0], SIX_V[2], SIX_V[4]
    for a2 in (e1, e2, -e0):
        for b, b2 in itertools.product(readouts, repeat=2):
            yield e0, a2, b, b2


def scan_hex(deltas, readouts=None):
    quads = list(quadruples_h(readouts))
    return {d: max(lp_max_chsh_h(a, a2, b, b2, d) for a, a2, b, b2 in quads) for d in deltas}


def main():
    print("closed timelines n = 0..7  hex  :", [closed_count_hex(n) for n in range(8)])
    print("                           cubic:", [closed_count_cubic(n) for n in range(8)])
    print("minimal closed words (ordered), hex   n=2,3,4:", [minimal_closed_words(n, HINT) for n in (2, 3, 4)])
    print("                                cubic n=2,3,4:", [minimal_closed_words(n, CUBIC_INT) for n in (2, 3, 4)])
    a, b = hex_measure_a(["+A"] + ["+B"] * 7, ["-A"] + ["+B"] * 7, 5, "+C")
    print("copy contract, A measured +C at step 5:\n  A =", a, "\n  B =", b)
    print("max CHSH on the six directions: toy %.3f, singlet %.3f (2 sqrt2 = %.3f)" % (chsh_max_six(toy_E), chsh_max_six(singlet_E), 2 * math.sqrt(2)))
    print("delta  max CHSH under no signalling (b among the 12 directions)")
    for d, s in scan_hex([0.0, 0.05, 0.1, 0.2, 0.3, 0.5, 5 / 6]).items():
        print("%.3f   %.4f" % (d, s))


if __name__ == "__main__":
    main()
