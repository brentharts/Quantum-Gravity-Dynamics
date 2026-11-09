#!/usr/bin/env python3
"""
contract_options.py -- three ways a "consistency contract" between an entangled pair of six-direction timelines could work, compared exactly.

Setting.  A measurement of A forces one step of A to s*alpha: alpha is one of the six directions (A's setting) and the outcome s = +-1 is fair.  A reads
out along alpha only.  B is read out along any of 26 directions beta (the 26 non-zero vectors of {-1,0,1}^3) by the rule outcome = sgn(beta . F_B), F_B being
B's net flux (0 when orthogonal).  A's other n-1 steps are uniform random.  n = 6 here.  The readout rule is an assumption of ours; the toy has none.

  G  global conservation: the pair's joint flux is kept at zero, and B repairs whatever A's forced step changed, so F_B = -(G + s alpha) with G the sum of A's other steps.
  L  local conservation: each particle keeps its own flux; B is never touched.
  C  copy (the toy as written): B is overwritten with -s*alpha and read at its last step, so E = -sgn(alpha . beta).
  R  revelation: nothing is forced; each side reads a pre-existing origin axis v (uniform over the six directions), A gets sgn(alpha.v) and B gets sgn(beta.(-v)).

Results (all Exact: the distributions are computed by convolution, not sampled)
  C  signals completely (B's history is -s*alpha) and reaches CHSH = 4.
  G  B's flux distribution depends on A's setting (total variation distance between alpha = X and alpha = Y reported), so the history signals.
  L  B's distribution does not depend on alpha, and the correlation is exactly 0, because A's outcome is independent of anything B holds.
  R  no signalling, and CHSH is bounded by 2 (Bell); the exact maximum over A in six directions and B in 26 is reported.
Together these are the four policies considered here: they signal (C, G), do not correlate (L), or are local hidden-variable models (R). The CHSH value of G and the
value 4 of C depend on how B is read (flux against last step), so only the signalling statements and the bounds of L and R are robust.
This is a statement about these policies, not a no-go for retrocausal models in general: a model whose history ensemble depends on both settings
(settings as future boundary conditions) is outside this class and is not built here.

    python3 directional/contract_options.py
"""
import itertools
from collections import Counter

import numpy as np

from axis_timelines import VEC

SIX = [np.array(v, dtype=int) for v in VEC.values()]
TWENTY_SIX = [np.array(v, dtype=int) for v in itertools.product((-1, 0, 1), repeat=3) if any(v)]


def flux_distribution(n):
    """Exact distribution {flux vector: probability} of the sum of n uniform six-direction steps."""
    dist = {(0, 0, 0): 1.0}
    for _ in range(n):
        nxt = Counter()
        for f, p in dist.items():
            for v in SIX:
                nxt[(f[0] + v[0], f[1] + v[1], f[2] + v[2])] += p / 6.0
        dist = dict(nxt)
    return dist


def b_flux_distribution(alpha, n):
    """G: distribution of F_B = -(G + s alpha) over fair s and G = sum of A's other n-1 steps."""
    g = flux_distribution(n - 1)
    out = Counter()
    for s in (+1, -1):
        for f, p in g.items():
            out[tuple(-(f[i] + s * alpha[i]) for i in range(3))] += p / 2.0
    return dict(out)


def total_variation(p, q):
    keys = set(p) | set(q)
    return 0.5 * sum(abs(p.get(k, 0.0) - q.get(k, 0.0)) for k in keys)


def correlation_G(alpha, beta, n=6):
    g = flux_distribution(n - 1)
    tot = 0.0
    for s in (+1, -1):
        for f, p in g.items():
            fb = [-(f[i] + s * alpha[i]) for i in range(3)]
            tot += 0.5 * p * s * np.sign(np.dot(beta, fb))
    return tot


def correlation_L(alpha, beta, n=6):
    return 0.0


def correlation_C(alpha, beta, n=6):
    return -float(np.sign(np.dot(alpha, beta)))


def correlation_R(alpha, beta, n=6):
    return -sum(np.sign(np.dot(alpha, v)) * np.sign(np.dot(beta, v)) for v in SIX) / 6.0


def chsh_max(E, n=6):
    """Maximum |E(a,b) - E(a,b') + E(a',b) + E(a',b')| over a, a' in the six directions and b, b' in the 26."""
    M = np.array([[E(a, b, n) for b in TWENTY_SIX] for a in SIX])
    best = 0.0
    for i, j in itertools.product(range(len(SIX)), repeat=2):
        x, y = M[i], M[j]
        s = np.abs((x[:, None] - x[None, :]) + (y[:, None] + y[None, :]))
        best = max(best, s.max())
    return best


def main(n=6):
    ex, ey = SIX[0], SIX[2]
    print("G: total variation of B's flux distribution, A's setting X against Y (n = %d): %.4f" %
          (n, total_variation(b_flux_distribution(ex, n), b_flux_distribution(ey, n))))
    print("G: same for X against -X: %.4f" % total_variation(b_flux_distribution(ex, n), b_flux_distribution(-ex, n)))
    for name, E in (("C", correlation_C), ("G", correlation_G), ("L", correlation_L), ("R", correlation_R)):
        print("%s: max CHSH, A in six directions, B in 26: %.4f" % (name, chsh_max(E, n)))


if __name__ == "__main__":
    main()
