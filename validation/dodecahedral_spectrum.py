#!/usr/bin/env python3
"""
dodecahedral_spectrum.py -- Laplace spectrum of the Poincare dodecahedral space S^3 / I* (exploratory).
=====================================================================================================

I* is the binary icosahedral group (order 120), built here from its unit quaternions.  Harmonic polynomials of degree k
on R^4 carry (k+1) copies of the spin-k/2 representation under left multiplication, so the number of I*-invariant modes is
(k+1) N_k, with N_k the multiplicity of the trivial representation in Sym^k(C^2)|I*, obtained by character averaging,
chi_k(g) = sin((k+1) theta) / sin(theta).  N_k reproduces Klein's Molien series (1 + t^30) / ((1 - t^12)(1 - t^20)).
Eigenvalue of -Laplacian on unit S^3: k (k + 2).  Full S^3: lowest nonzero k = 1 (eigenvalue 3); here k = 12 (168).
No C_l is computed, and QGD makes no statement about global topology.
"""
import itertools

import numpy as np

PHI = (1 + 5 ** 0.5) / 2


def _even_permutations(n=4):
    for p in itertools.permutations(range(n)):
        inv = sum(1 for i in range(n) for j in range(i + 1, n) if p[i] > p[j])
        if inv % 2 == 0:
            yield p


def binary_icosahedral_group():
    """The 120 unit quaternions (w, x, y, z) of I*."""
    els = []
    for i in range(4):                                   # +-1, +-i, +-j, +-k
        for s in (1, -1):
            q = [0.0] * 4; q[i] = float(s); els.append(q)
    for signs in itertools.product((0.5, -0.5), repeat=4):  # (+-1 +-i +-j +-k)/2
        els.append(list(signs))
    base = (PHI / 2, 0.5, 1 / (2 * PHI), 0.0)            # even permutations of (phi, 1, 1/phi, 0)/2, all signs
    for p in _even_permutations():
        v = [base[p[i]] for i in range(4)]
        nz = [i for i in range(4) if v[i] != 0.0]
        for signs in itertools.product((1, -1), repeat=3):
            q = list(v)
            for i, s in zip(nz, signs):
                q[i] = s * v[i]
            els.append(q)
    return np.array(els)


def qmul(a, b):
    w1, x1, y1, z1 = a; w2, x2, y2, z2 = b
    return np.array([w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2, w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
                     w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2, w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2])


def is_closed(G, tol=1e-9):
    """True if the set is closed under multiplication (all 120^2 products are members)."""
    G = np.asarray(G)
    for a in G:
        prods = np.array([qmul(a, b) for b in G])
        d = np.abs(prods[:, None, :] - G[None, :, :]).sum(axis=2).min(axis=1)
        if d.max() > tol:
            return False
    return True


def invariants(G, kmax):
    """N_k = (1/|G|) sum_g chi_k(g) for k = 0..kmax, by character averaging."""
    G = np.asarray(G)
    theta = np.arccos(np.clip(G[:, 0], -1, 1))
    N = np.zeros(kmax + 1)
    for k in range(kmax + 1):
        s = np.sin(theta)
        chi = np.where(s > 1e-12, np.sin((k + 1) * theta) / np.where(s > 1e-12, s, 1.0),
                       np.where(np.cos(theta) > 0, k + 1.0, (-1.0) ** k * (k + 1.0)))
        N[k] = chi.sum() / len(G)
    return N


def molien_coefficients(kmax):
    """Coefficients of (1 + t^30) / ((1 - t^12)(1 - t^20)) up to t^kmax."""
    c = np.zeros(kmax + 1, dtype=int)
    for a in range(0, kmax + 1, 12):
        for b in range(0, kmax + 1 - a, 20):
            c[a + b] += 1
            if a + b + 30 <= kmax:
                c[a + b + 30] += 1
    return c


def spectrum(G, kmax):
    """List of (k, eigenvalue k(k+2), multiplicity (k+1) N_k) for the k <= kmax with N_k > 0."""
    N = np.rint(invariants(G, kmax)).astype(int)
    return [(k, k * (k + 2), (k + 1) * int(N[k])) for k in range(kmax + 1) if N[k] > 0]


if __name__ == "__main__":
    G = binary_icosahedral_group()
    print("order", len(G), "closed", is_closed(G))
    sp_ = spectrum(G, 60)
    for k, ev, m in sp_[:8]:
        print(f"k={k:3d}  eigenvalue={ev:5d}  multiplicity={m}")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.vlines([s[1] for s in sp_ if s[0] > 0], 0, [s[2] for s in sp_ if s[0] > 0], label="S^3 / I*")
    ax.axvline(3, color="r", ls="--", label="S^3 lowest (k=1)")
    ax.set_xlabel("eigenvalue k(k+2)"); ax.set_ylabel("multiplicity"); ax.legend()
    fig.tight_layout(); fig.savefig(__file__.replace(".py", ".png"), dpi=120)
