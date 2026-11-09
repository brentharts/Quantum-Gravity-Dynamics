#!/usr/bin/env python3
"""
hex_contract.py -- repairing a hexagonal timeline after a shock, with and without a reward for closed 3-cycles.

The consistency contract of the toy (axis_timelines.py) edits steps until the net flux is zero.  With the six hexagonal directions (hex_directions.py) a
timeline can also close in three steps (e_0 + e_2 - e_1 = 0).  This module asks what that does to the repair of a shocked timeline, exactly for the
minimum number of edits, and numerically for a greedy solver whose cost can reward locally closed triples.

Exact C1.  Minimum repair.  Take a closed timeline of n steps, force one step to a different direction (the shock), and ask for the fewest OTHER steps to change
           so that the timeline is closed again (exhaustive over all closed timelines, all shock positions and all new directions).
             cubic, n = 6 :  always exactly 1 other step (1860 closed timelines x 6 x 5 shocks);
             hexagonal, n = 5 :  1 other step in 7560 of 9000 shocks, 2 in 1440;   n = 6 : 1 in 57960 of 61200, 2 in 3240;   n = 7 : 1 in 335160 of 352800, 2 in 17640; never more than 2.
           So hexagonal closure is not cheaper to repair: the 3-cycles give more closed timelines but a single shock can need two coordinated edits.
Exact C2.  Odd lengths.  A cubic timeline of odd length has no closed completion at all (no closed words exist), so no contract can repair it to zero flux; the
           hexagonal set has closed timelines of every length >= 2 (n = 3: 12, n = 5: 360, n = 7: 10080).
Numeric C3.  Greedy solver (one sweep = try every non-shocked step, take the direction with the lowest cost, stop when nothing improves), seeded.  Cost =
           |flux|^2 - mu * T, T = number of consecutive triples that sum to zero.  Over 800 seeded trials per row (closed timeline, one random shock):
                           n = 6                          n = 8
                        closes  edits  reach  T            closes  edits  reach  T
             cubic         1.00    1.13   2.4   0.0           1.00    1.14   3.2   0.0
             hex  mu=0     0.92    1.22   2.4   0.2           0.96    1.28   3.2   0.5
             hex  mu=.5    0.92    1.49   2.6   0.8           0.98    1.55   3.5   0.9
             hex  mu=1.5   0.73    2.41   3.2   2.2           0.53    3.42   4.7   3.2
           (closes = fraction ending with zero flux; edits = distinct steps changed; reach = largest distance from the shock of a changed step; T = triples at the end.)
           Reading: rewarding 3-cycles makes the final history built from locally closed triples (T rises) and pays for it: the retro-wave edits more steps and reaches
           further, and with a strong reward (mu = 1.5) the solver trades global closure for local closure and often stops with non-zero flux.  The zero-reward hexagonal
           solver also gets stuck (0.92-0.96), because its single-step moves cannot do the two coordinated edits that C1 says are sometimes needed.
What this is not: a derivation of the reward from E8 or from the tile (it is a choice of cost), a statement about dynamics, or about the sign of the arrow of time.

    python3 directional/hex_contract.py
"""
import itertools
import random
from collections import Counter

import numpy as np

HEX = [(1, 0), (-1, 0), (0, 1), (0, -1), (-1, 1), (1, -1)]
CUB = [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)]
_C = {}


def closed_words(vec, n):
    key = (tuple(vec), n)
    if key not in _C:
        V = np.array(vec)
        out = [w for w in itertools.product(range(6), repeat=n) if not V[list(w)].sum(0).any()]
        _C[key] = np.array(out, dtype=np.int8).reshape(len(out), n)
    return _C[key]


def repair_distribution(vec, n):
    """Counter {extra edits: number of shocks} over all closed timelines, shock positions and new directions; key -1 = no closed completion."""
    C = closed_words(vec, n)
    cnt = Counter()
    if len(C) == 0:
        return cnt
    for w in C:
        for t in range(n):
            for s in range(6):
                if s == w[t]:
                    continue
                sel = C[C[:, t] == s]
                cnt[-1 if len(sel) == 0 else int((sel != w).sum(1).min()) - 1] += 1
    return cnt


# ---------------------------------------------------------------- greedy solver
def _norm2(f, kind):
    return f[0] * f[0] + f[0] * f[1] + f[1] * f[1] if kind == "hex" else sum(x * x for x in f)


def _flux(w, V):
    return tuple(sum(V[a][i] for a in w) for i in range(len(V[0])))


def triples(w, V):
    return sum(1 for i in range(len(w) - 2) if not any(V[w[i]][k] + V[w[i + 1]][k] + V[w[i + 2]][k] for k in range(len(V[0]))))


def cost(w, V, kind, mu):
    return _norm2(_flux(w, V), kind) - mu * triples(w, V)


def solve(w, shock, V, kind, mu, rng, max_sweeps=30):
    w = list(w)
    n = len(w)
    edits = set()
    for _ in range(max_sweeps):
        changed = False
        order = [t for t in range(n) if t != shock]
        rng.shuffle(order)
        for t in order:
            best, pick = cost(w, V, kind, mu), w[t]
            for a in range(6):
                if a == w[t]:
                    continue
                old = w[t]
                w[t] = a
                c = cost(w, V, kind, mu)
                w[t] = old
                if c < best - 1e-12:
                    best, pick = c, a
            if pick != w[t]:
                w[t] = pick
                edits.add(t)
                changed = True
        if not changed:
            break
    return w, edits


def run(kind, n, mu, trials=800, seed=3):
    V = HEX if kind == "hex" else CUB
    rng = random.Random(seed)
    ok, ne, far, tri = 0, [], [], []
    for _ in range(trials):
        while True:
            w = [rng.randrange(6) for _ in range(n)]
            if not any(_flux(w, V)):
                break
        s = rng.randrange(n)
        w[s] = rng.choice([a for a in range(6) if a != w[s]])
        w2, ed = solve(w, s, V, kind, mu, rng)
        ok += not any(_flux(w2, V))
        ne.append(len(ed))
        far.append(max([abs(t - s) for t in ed], default=0))
        tri.append(triples(w2, V))
    return dict(closes=ok / trials, edits=float(np.mean(ne)), reach=float(np.mean(far)), T=float(np.mean(tri)))


def main():
    for n in (5, 6):
        print("hex   n=%d" % n, dict(repair_distribution(HEX, n)))
    print("cubic n=6", dict(repair_distribution(CUB, 6)))
    print("cubic n=5 (no closed words):", dict(repair_distribution(CUB, 5)))
    for n in (6, 8):
        for kind, mu in [("cub", 0), ("hex", 0), ("hex", 0.5), ("hex", 1.5)]:
            print(n, kind, mu, {k: round(v, 3) for k, v in run(kind, n, mu).items()})


if __name__ == "__main__":
    main()
