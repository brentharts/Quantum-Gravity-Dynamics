#!/usr/bin/env python3
"""
update_rate_clock.py -- clock rate and redshift from the toy's own update rule.

In GravitationalTimelineNetwork a node of density rho takes part in a relaxation sweep with probability 1/rho (it is skipped when random() > 1/rho).  Read that as a
clock: a node ticks when it updates, and the sweep count is the coordinate time.  Then

Exact 1.  The expected number of ticks in T sweeps is T/rho and the variance is T (1/rho)(1 - 1/rho), so the clock rate is 1/rho.  This follows from the update
          rule; it is not an extra assumption beyond that rule.
Exact 2.  Station 1 emits a pulse every K of its own ticks.  The pulse takes a fixed number of sweeps to arrive (static setting), so station 2 counts, in expectation,
          K rho_1 / rho_2 ticks between receptions: the ratio of received to emitted intervals is clock_2 / clock_1, the gravitational redshift between the stations.
Exact 3.  With rho = 1 / sqrt(1 - r_s / r) (the QGD-matching density of clock_rate_vs_sigma.py) this ratio is sqrt((1 - r_s/r_2) / (1 - r_s/r_1)), the Schwarzschild redshift.
Numeric.  Seeded simulation of both statements, with the sampling error.  Nothing here derives rho: the density is still an input, and what is derived is that an
          update probability 1/rho gives a clock 1/rho and a redshift equal to the density ratio.

    python3 directional/update_rate_clock.py
"""
import math

import numpy as np


def rho_qgd(r_over_rs):
    return 1.0 / math.sqrt(1.0 - 1.0 / r_over_rs)


def simulate_ticks(rho, sweeps, rng):
    """Number of ticks in `sweeps` sweeps for a node updating with probability 1/rho per sweep."""
    return int(np.count_nonzero(rng.random(sweeps) < 1.0 / rho))


def tick_statistics(rho, sweeps):
    p = 1.0 / rho
    return sweeps * p, math.sqrt(sweeps * p * (1.0 - p))


def received_ticks(rho1, rho2, k_emitted, rng):
    """Station 1 emits after k_emitted of its own ticks: the waiting time in sweeps is a sum of geometric variables; station 2 then counts its own ticks in that time."""
    waiting_sweeps = int(rng.negative_binomial(k_emitted, 1.0 / rho1) + k_emitted)      # sweeps until the k-th tick
    return simulate_ticks(rho2, waiting_sweeps, rng)


def redshift_prediction(r1_over_rs, r2_over_rs):
    return math.sqrt((1.0 - 1.0 / r2_over_rs) / (1.0 - 1.0 / r1_over_rs))


def main(seed=11):
    rng = np.random.default_rng(seed)
    r1, r2 = 3.0, 10.0
    rho1, rho2 = rho_qgd(r1), rho_qgd(r2)
    T = 2_000_000
    n1, n2 = simulate_ticks(rho1, T, rng), simulate_ticks(rho2, T, rng)
    print("rho(r/rs = %.0f) = %.4f, rho(%.0f) = %.4f" % (r1, rho1, r2, rho2))
    for name, rho, n in (("station 1", rho1, n1), ("station 2", rho2, n2)):
        mean, sd = tick_statistics(rho, T)
        print("%s: %d ticks in %d sweeps, expected %.0f +- %.0f (%.1f sigma)" % (name, n, T, mean, sd, (n - mean) / sd))
    K = 200
    got = np.mean([received_ticks(rho1, rho2, K, rng) for _ in range(4000)]) / K
    print("ticks of station 2 per tick of station 1: simulated %.4f, rho1/rho2 = %.4f, Schwarzschild sqrt((1-rs/r2)/(1-rs/r1)) = %.4f"
          % (got, rho1 / rho2, redshift_prediction(r1, r2)))


if __name__ == "__main__":
    main()
