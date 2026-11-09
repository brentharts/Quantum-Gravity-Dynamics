#!/usr/bin/env python3
"""
axis_timelines.py -- the six-direction timeline model: exact counts, and what the "consistency contract" does.

The toy models in spectre-monotile-py/research/toymodel_v1 (AxisNetwork.py and relatives) describe a particle history as a list
of steps, each one of the six signed axis directions +X, -X, +Y, -Y, +Z, -Z.  A "retrocausal consistency contract" then edits
steps before and after a shocked step until the vector sum of the timeline (its net flux) is small.  This module restates that
model with a seeded random generator, so every number can be reproduced, and checks what can be checked exactly.

Labels (as in the rest of this repository):  Exact = enumeration or identity, tested;  Numeric = computed with a stated
sample and seed;  Hypothesis = not shown.  Nothing here says anything about spacetime or about quantum mechanics.

Results
-------
Exact 1.  The number of timelines of n steps with zero net flux (closed walks on Z^3) is 1, 0, 6, 0, 90, 0, 1860, 0, 44730 for
          n = 0..8.  The generating-function count here is checked against brute force for n <= 8 (6^8 = 1 679 616 timelines).
Exact 2.  |flux|^2 has the same parity as n.  So for an even number of steps the contract's stopping rule |flux|^2 <= 3 means
          |flux|^2 in {0, 2}, and for an odd number |flux|^2 in {1, 3}.  (Brute force for n <= 8.)
Exact 3.  The edit that the repository's solver intends, "replace step t by the step that makes the new flux smallest", changes the
          flux by  new - old.  The solver as written (AxisNetwork._resolve_consistency) scores the candidate by
          |flux + old + new|^2, i.e. with old added instead of removed.  `solver_cost_matches_flux_after_edit` tests the difference.
Exact 4.  In the gravitational relaxation (GravitationalTimelineNetwork) the update copies step t-1 into step t, so with uniform density 1 (no
          node skipped) the forced event is overwritten in the first sweep and steps 1..n-2 all equal step 0 afterwards.
Numeric.  Over seeded random timelines the two versions (as written, and with the sign corrected) are compared by how often they
          reach |flux|^2 <= 3 within the iteration limit.

    python3 directional/axis_timelines.py          # prints the tables
"""
import itertools
import random
from collections import Counter

DIRS = ("+X", "-X", "+Y", "-Y", "+Z", "-Z")
VEC = {"+X": (1, 0, 0), "-X": (-1, 0, 0), "+Y": (0, 1, 0), "-Y": (0, -1, 0), "+Z": (0, 0, 1), "-Z": (0, 0, -1)}


def flux(timeline):
    """Net flux: the vector sum of the steps (the endpoint of the walk)."""
    return tuple(sum(VEC[a][i] for a in timeline) for i in range(3))


def flux2(timeline):
    return sum(x * x for x in flux(timeline))


def closed_count(n):
    """Exact number of timelines of n steps with zero net flux, by dynamic programming over positions."""
    pos = Counter({(0, 0, 0): 1})
    for _ in range(n):
        nxt = Counter()
        for p, c in pos.items():
            for a in DIRS:
                v = VEC[a]
                nxt[(p[0] + v[0], p[1] + v[1], p[2] + v[2])] += c
        pos = nxt
    return pos[(0, 0, 0)]


def brute_flux_counts(n):
    """{|flux|^2: number of timelines} over all 6^n timelines, by enumeration."""
    out = Counter()
    for t in itertools.product(DIRS, repeat=n):
        out[flux2(t)] += 1
    return out


def solver_cost_matches_flux_after_edit(timeline, t, ax):
    """(cost used by the solver as written, true |flux|^2 after replacing step t by ax) for a given edit."""
    f = flux(timeline)
    old = VEC[timeline[t]]
    new = VEC[ax]
    as_written = sum((old[i] + new[i] + f[i]) ** 2 for i in range(3))
    edited = list(timeline)
    edited[t] = ax
    return as_written, flux2(edited)


def resolve(timeline, shock, corrected=False, max_iter=30):
    """The consistency contract of AxisNetwork._resolve_consistency, without printing.  Returns (timeline, iterations, converged).
    corrected=False reproduces the cost  |old + new + flux|^2  as written; corrected=True uses |flux - old + new|^2."""
    tl = list(timeline)
    for it in range(max_iter):
        f = flux(tl)
        if sum(x * x for x in f) <= 3:
            return tl, it, True
        spread = it // 3 + 1
        for off in range(-spread, spread + 1):
            tgt = shock + off
            if 0 <= tgt < len(tl) and tgt != shock:
                old = VEC[tl[tgt]]
                sign = -1 if corrected else 1
                tl[tgt] = min(DIRS, key=lambda ax: sum((sign * old[i] + VEC[ax][i] + f[i]) ** 2 for i in range(3)))
    return tl, max_iter, sum(x * x for x in flux(tl)) <= 3


def gravity_relax(timeline, density, shock, new_axis, rng, max_iter=25):
    """GravitationalTimelineNetwork.apply_gravitational_relaxation without printing: force step `shock`, then sweep t = 1..n-2, each node
    skipped with probability 1 - 1/density[t], otherwise set to the direction closest to step t-1 (which is step t-1 itself)."""
    tl = list(timeline)
    tl[shock] = new_axis
    for _ in range(max_iter):
        changed = False
        for t in range(1, len(tl) - 1):
            if rng.random() > 1.0 / density[t]:
                continue
            ref = VEC[tl[t - 1]]
            best = min(DIRS, key=lambda ax: sum((VEC[ax][i] - ref[i]) ** 2 for i in range(3)))
            if tl[t] != best:
                tl[t] = best
                changed = True
        if not changed:
            break
    return tl


def convergence_table(n=12, trials=2000, seed=1):
    """Numeric: fraction of random timelines (after forcing one random step to a random direction) that the solver brings to |flux|^2 <= 3."""
    rng = random.Random(seed)
    rows = {}
    for corrected in (False, True):
        ok = 0
        its = []
        r = random.Random(seed)
        for _ in range(trials):
            tl = [r.choice(DIRS) for _ in range(n)]
            k = r.randrange(n)
            tl[k] = r.choice(DIRS)
            out, it, conv = resolve(tl, k, corrected)
            ok += conv
            its.append(it)
        rows["corrected" if corrected else "as written"] = (ok / trials, sum(its) / trials)
    return rows


def main():
    print("closed timelines (zero net flux) of n steps; dynamic programming against brute force")
    for n in range(0, 9):
        dp = closed_count(n)
        bf = brute_flux_counts(n)[0] if n <= 8 else None
        print("  n = %d  %6d  brute force %6d" % (n, dp, bf))
    print("\n|flux|^2 values that occur, n = 8 (parity of n):", sorted(brute_flux_counts(8)))
    print("\nconsistency contract, n = 12, 2000 seeded timelines, one random step forced:")
    for k, (frac, mean_it) in convergence_table().items():
        print("  %-11s reaches |flux|^2 <= 3 in %.1f%% of runs, mean iterations %.1f" % (k, 100 * frac, mean_it))


if __name__ == "__main__":
    main()
