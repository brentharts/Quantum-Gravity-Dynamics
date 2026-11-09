#!/usr/bin/env python3
"""
hex_e8.py -- the six hexagonal directions as roots of E8, and the E8 read-outs of the Bell test.

hex_directions.py uses six directions 60 degrees apart (the red A2 roots of the Spectre lift) but only as a direction set.  Here they are put into E8 by the
embedding of spectre_in_e8 (spectre-monotile-py, spectre repository): E8 as Hurwitz pairs (a, b) with a - b in (1+i)H, fibre rotation R_w for
w = (-1+i+j+k)/2 (order 3, no fixed vector), r1 = (1+i, 0), and the red lift (p, q) goes to  p r1 + q (r1 + R_w r1).  (Written here from the description, not
imported: this repository does not depend on the others.)

Results (Exact; tests in tests/test_directional_hex_e8.py)
Exact E1.  The six unit hexagonal vectors map to six roots of E8 (norm^2 2) with <embed x, embed y> = 2 <x, y> for the hexagonal inner product, i.e. the
           map is an isometry of A2 onto an A2 root subsystem of E8, so "sign of the hexagonal inner product" = "sign of the E8 inner product".
Exact E2.  The roots of E8 orthogonal to this A2 number 72 (an E6), and the remaining 240 - 6 - 72 = 162 roots have inner products +-1 with some of the six.
Exact E3.  The E8-intrinsic read-out set is the 240 roots.  Their sign patterns on the six hexagonal roots (sgn <b, v>, v in the six) take a small number of
           values; roots with a zero in the pattern exist, so, unlike the planar 12 directions of hex_directions.py, the intrinsic read-outs have ternary outcomes.
           The red and green roots are orthogonal in E8 (tested in spectre_e8.py of the spectre repository, not here), so a green root reads 0 on every red
           hidden state: the planar picture (green 30 degrees from red, non-zero inner products) is NOT the E8 picture.
Numeric E4.  Price of retro-dependence with E8 read-outs (hidden state = one of the six red roots, a among the six, b and b' among the 13 distinct sign
           patterns of the 240 roots), by the LP of retro_boundary_model.py: the largest CHSH under no signalling is min(4, 2 + 12 delta), the same as for the planar
           read-outs of hex_directions.py (checked at delta = 0, 0.05, 0.1, 0.2).  The ternary patterns do not raise it.

    python3 directional/hex_e8.py
"""
import itertools

import numpy as np
from scipy.optimize import linprog

_C = {}
OMEGA = np.array([-0.5, 0.5, 0.5, 0.5])


def qmul(p, q):
    a1, b1, c1, d1 = p
    a2, b2, c2, d2 = q
    return np.array([a1 * a2 - b1 * b2 - c1 * c2 - d1 * d2, a1 * b2 + b1 * a2 + c1 * d2 - d1 * c2,
                     a1 * c2 - b1 * d2 + c1 * a2 + d1 * b2, a1 * d2 + b1 * c2 - c1 * b2 + d1 * a2])


def _in_hurwitz(y):
    return bool(np.allclose(y, np.round(y), atol=1e-9) or np.allclose(y - 0.5, np.round(y - 0.5), atol=1e-9))


def roots():
    if "R" not in _C:
        U = [np.eye(4)[i] * s for i in range(4) for s in (1.0, -1.0)] + [np.array(t) for t in itertools.product((0.5, -0.5), repeat=4)]
        omi = np.array([1.0, -1.0, 0, 0])
        R = [np.r_[a, b] for a in U for b in U if _in_hurwitz(qmul(omi, a - b) / 2.0)]
        n2 = [qmul(np.array([1.0, 1, 0, 0]), u) for u in U]
        R += [np.r_[np.zeros(4), b] for b in n2] + [np.r_[a, np.zeros(4)] for a in n2]
        _C["R"] = np.array(R)
    return _C["R"]


def right_matrix(l):
    M = np.zeros((8, 8))
    for k in range(4):
        e = np.zeros(4)
        e[k] = 1
        c = qmul(e, l)
        M[k, :4] = c
        M[4 + k, 4:] = c
    return M


R1 = np.r_[1.0, 1, 0, 0, np.zeros(4)]


def embed_red(p, q):
    """E8 vector of the red lift coordinates (p, q): p r1 + q (r1 + R_w r1)."""
    Rw = _C.setdefault("Rw", right_matrix(OMEGA))
    return p * R1 + q * (R1 + R1 @ Rw)


# integer coordinates (p, q) of the six hexagonal unit vectors e^{i pi k/3}, in the order +A, -A, +B, -B, +C, -C of hex_directions.py:
# +A = 1 = (1,0); +B = e^{i60} = w6 = (0,1); +C = e^{i120} = w6^2 = w6 - 1 = (-1,1)
HEX_PQ = {"+A": (1, 0), "-A": (-1, 0), "+B": (0, 1), "-B": (0, -1), "+C": (-1, 1), "-C": (1, -1)}


def hex_roots():
    return np.array([embed_red(*pq) for pq in HEX_PQ.values()])


def orthogonal_roots():
    A = hex_roots()
    return [r for r in roots() if all(abs(r @ a) < 1e-9 for a in A)]


def sign_pattern(b, A=None):
    A = hex_roots() if A is None else A
    return tuple(int(np.sign(round(float(b @ a)))) for a in A)


def distinct_readouts():
    """Distinct sign patterns of the 240 roots on the six hexagonal roots -> {pattern: representative root}."""
    out = {}
    for r in roots():
        out.setdefault(sign_pattern(r), r)
    return out


# ---------------------------------------------------------------- LP with E8 read-outs
def _lp(a_idx, a2_idx, b_vec, b2_vec, delta):
    A = hex_roots()
    n = 6

    def oa(i, v):
        return int(np.sign(round(float(A[i] @ A[v]))))

    def ob(b, v):
        return int(np.sign(round(float(b @ -A[v]))))
    pairs = [(("a", a_idx), ("b", b_vec)), (("a", a_idx), ("b", b2_vec)), (("a", a2_idx), ("b", b_vec)), (("a", a2_idx), ("b", b2_vec))]
    signs = [+1, -1, +1, +1]

    def corr(pa, pb):
        return np.array([oa(pa[1], v) * ob(pb[1], v) for v in range(n)], dtype=float)
    best = 0.0
    for overall in (+1, -1):
        c = np.zeros(4 * n)
        for k, ((pa, pb), sg) in enumerate(zip(pairs, signs)):
            c[k * n:(k + 1) * n] = -overall * sg * corr(pa, pb)
        A_eq, b_eq = [], []
        for k in range(4):
            row = np.zeros(4 * n)
            row[k * n:(k + 1) * n] = 1
            A_eq.append(row)
            b_eq.append(1.0)

        def add_equal(k1, k2, f):
            for s in (-1, 0, 1):
                row = np.zeros(4 * n)
                for v in range(n):
                    if f(v) == s:
                        row[k1 * n + v] += 1
                        row[k2 * n + v] -= 1
                A_eq.append(row)
                b_eq.append(0.0)
        add_equal(0, 1, lambda v: oa(a_idx, v))
        add_equal(2, 3, lambda v: oa(a2_idx, v))
        add_equal(0, 2, lambda v: ob(b_vec, v))
        add_equal(1, 3, lambda v: ob(b2_vec, v))
        bounds = [(max(0.0, 1 / 6 - delta), min(1.0, 1 / 6 + delta))] * (4 * n)
        res = linprog(c, A_eq=np.array(A_eq), b_eq=np.array(b_eq), bounds=bounds, method="highs")
        if res.status == 0:
            best = max(best, -res.fun)
    return best


def scan_e8(deltas):
    """{delta: largest CHSH}: a = hexagonal root 0, a' in {60, 120, 180 degrees}, b and b' over the distinct sign patterns of the 240 roots."""
    reps = list(distinct_readouts().values())
    quads = [(0, a2, b, b2) for a2 in (2, 4, 1) for b in reps for b2 in reps]
    return {d: max(_lp(*q, d) for q in quads) for d in deltas}


def main():
    A = hex_roots()
    G = A @ A.T
    print("Gram of the six embedded roots (E8 inner products):\n", np.round(G).astype(int))
    print("roots orthogonal to the A2:", len(orthogonal_roots()))
    pats = distinct_readouts()
    print("distinct sign patterns of the 240 roots on the six:", len(pats), "; with a zero entry:", sum(0 in p for p in pats))
    print("delta  max CHSH (E8 read-outs)")
    for d, s in scan_e8([0.0, 0.05, 0.1, 0.2]).items():
        print("%.3f   %.4f" % (d, s))


if __name__ == "__main__":
    main()
