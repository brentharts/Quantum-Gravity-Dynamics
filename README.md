# Quantum Gravitational Dynamics (QGD)

> **Thesis (Romeo Matshaba):** spacetime geometry is not fundamental. It is built algebraically from a dimensionless phase field $\sigma_\mu$ obtained from the semiclassical (WKB) limit of the Dirac spinor, and General Relativity is recovered as a consequence.

[![Theory](https://img.shields.io/badge/Theory-QGD-blue)](docs/complete_paper.tex)
[![Solutions](https://img.shields.io/badge/GR_Solutions-Algebraic-green)](comparison/EFE_solutions_from_QGD_perspective.py)
[![Dark Matter](https://img.shields.io/badge/Dark_Matter-No_Per--Galaxy_Fit-orange)](validation/)
[![License](https://img.shields.io/badge/License-MIT-lightgrey)](LICENSE)

This repository holds the theory, the code that implements it, and an independent test suite that checks it. Every
claim below carries a status taken from [`docs/claims_ledger.md`](docs/claims_ledger.md):
**verified** (re-derived independently), **recorded** (true with a stated limit), **fails** (false as stated),
**untested** (needs input not in this repository).

---

## Contents

1. [The idea](#the-idea)
2. [Master equations](#master-equations)
3. [Recovering GR solutions](#recovering-gr-solutions)
4. [Dark matter](#dark-matter)
5. [Claims and status](#claims-and-status)
6. [What the tests found](#what-the-tests-found)
7. [Open decisions for the author](#open-decisions-for-the-author)
8. [Repository structure](#repository-structure)
9. [Quick start](#quick-start)
10. [Citation](#citation)

---

## The idea

GR takes the metric $g_{\mu\nu}$ as fundamental. QGD takes the Dirac spinor $\psi = R\,e^{iS/\hbar}$, defines the
dimensionless phase gradient

$$\sigma_\mu \equiv \frac{1}{c}\,\partial_\mu S ,$$

and constructs the metric from it:

$$g_{\mu\nu} = T^\alpha_\mu T^\beta_\nu\Big(M_{\alpha\beta}\circ\big[\eta_{\alpha\beta} - \sum_{a}\varepsilon_a\,\sigma^{(a)}_\alpha\sigma^{(a)}_\beta - \kappa\ell_Q^2\,\partial_\alpha\sigma^\gamma\partial_\beta\sigma_\gamma\big]\Big)$$

with $\varepsilon_a=\pm1$ the source signature (attractive or repulsive), $\ell_Q$ the quantum gravitational length and
$\kappa\approx 2$ the stiffness. (The sign convention of this formula versus the code is an open item, see
[decisions](#open-decisions-for-the-author).)

```
Dirac spinor  psi = R exp(iS/hbar)
      | WKB limit
phase field   sigma_mu = (1/c) d_mu S
      | algebraic construction
metric        g = eta - sum eps sigma sigma - kappa l_Q^2 (d sigma)(d sigma)
      | variational principle
field eq.     box_g sigma_mu = Q_mu + G_mu + T_mu + kappa l_Q^2 box_g^2 sigma_mu
      | l_Q -> 0
Einstein      G_mn = (8 pi G / c^4) T_mn
```

The relation to GR is an equivalence under variable substitution, not a derivation from GR's field equations.

---

## Master equations

**Action**

$$S_\sigma = \int d^4x\sqrt{-g(\sigma)}\Big[-\frac{c^4}{16\pi G}R[g(\sigma)] + \tfrac12\nabla_\mu\sigma_\nu\nabla^\mu\sigma^\nu - \tfrac{\ell_Q^2}{2}\nabla_\alpha\nabla_\beta\sigma_\mu\nabla^\alpha\nabla^\beta\sigma^\mu\Big] + S_{\rm matter}$$

**Field equation** (fourth order)

$$\Box_g\sigma_\mu = Q_\mu + G_\mu + T_\mu + \kappa\ell_Q^2\,\Box_g^2\sigma_\mu + \mathcal O(\ell_Q^4)$$

$Q_\mu$ nonlinear self-interaction (the dark-matter sector), $G_\mu$ Kerr-Schild and radiative couplings,
$T_\mu=\tfrac12T^{\mu\nu}\sigma_\nu$ matter.

**Universal scalar.** All static and stationary metrics are written through one function, $g_{tt}=-(1-f)$:

$$f(r,\theta)=\underbrace{\tfrac{2GM}{c^2r}}_{\rm mass}-\underbrace{\tfrac{GQ^2}{c^4r^2}}_{\rm charge}+\underbrace{\tfrac{2Mr}{\Sigma}}_{\rm spin}+\underbrace{\tfrac{\Lambda r^2}{3}}_{\Lambda}+\underbrace{\tfrac{b(r)}{r}}_{\rm pressure}+\underbrace{H^2(t)r^2}_{\rm expansion}-\underbrace{\int\tfrac{P}{\rho c^2}dr}_{\rm EOS}+\underbrace{\kappa\tfrac{\hbar^2}{M^2c^2r^2}}_{\rm quantum}$$

In the static form $G^t{}_t$ is affine in $f$, so every term has a definite implied energy density
$\rho=(rf)'/(8\pi G r^2)$ (see [what the tests found](#what-the-tests-found)).

---

## Recovering GR solutions

The source-signature recipe: pick $\sigma^{(a)}$ and $\varepsilon_a$, apply the master metric.

| Spacetime | $\sigma$-field | $\varepsilon$ | Status |
|---|---|---|---|
| Schwarzschild | $\sqrt{2GM/c^2r}$ | $+1$ | verified (Kerr-Schild form, Ricci-flat) |
| Reissner-Nordström | $\sqrt{GQ^2/c^4r^2}$ | $-1$ | verified (exactly the Maxwell stress) |
| de Sitter / Nariai | $Hr$ | $+1$ | verified (double horizon at $r=1/\sqrt\Lambda$, area $4\pi/\Lambda$, $GM=c^3/3\sqrt3H$) |
| Anti-de Sitter | $\lvert H\rvert r$ | $-1$ | not separately tested |
| Kerr | $a\sin\theta\sqrt{2GM/c^2r}$ | $+1$ | untested; $g_{t\phi}$ from the cross term is reproduced, $g_{rr},g_{\theta\theta}$ are not supplied by $\sigma$ |

Frame dragging appears as interference: the Kerr $g_{t\phi}$ comes from the product of the mass and spin $\sigma$'s.
[`comparison/EFE_solutions_from_QGD_perspective.py`](comparison/EFE_solutions_from_QGD_perspective.py) builds ten of
these and prints `Passed 10/10 checks`.

**Limit of the recipe.** It is exact when all terms share one null direction (mass + charge + $\Lambda$ on one centre:
for a Kerr-Schild metric the mixed Ricci tensor is linear in $h$). It is **not** exact for separated sources: the
$\sigma_t=\sum_a\sqrt{2GM_a/c^2r_a}$ N-body metric has a nonzero Ricci tensor of order $(GM)^2$ from the cross term.
So the "exact N-body solution" is a statement about the $\sigma$-field sum, not a vacuum solution of Einstein's equations.

---

## Dark matter

The Taylor expansion of the phase factor gives factorial enhancement constants

$$\kappa_j=\sqrt{\frac{(2j-1)!}{2^{2j-2}}}\ \Rightarrow\quad \kappa = [1.00,\ 1.225,\ 2.74,\ 8.87,\ 37.7,\ 197,\ 1233, \ldots]$$

and a modified rotation curve $v^2=\frac{GM_{\rm bar}(<r)}{r}\big[1+\sum_{j=2}^{4}\kappa_j g_j(r)\big]$.

The engine is [`validation/dark_matter.py`](validation/dark_matter.py). It has no per-galaxy fit parameters but does
have shared empirical constants (`SIGMA_CRIT`, `ALPHA`, `BETA_0`, `LOG_M_TRIGGER`, two mass-to-light ratios, `G_CRIT`),
so it is not parameter-free.

| SPARC quantity | Author's value | Reproduced here | Note |
|---|---|---|---|
| pooled $R^2$ | 0.908 | 0.91 | verified |
| RMS | 8.3 km/s | about 26 km/s | fails under the choices in `reproduce_sparc.py` |
| $\chi^2_\nu$ | 1.2 | about 150 | fails, same choices |
| data set | 4,248 points, 175 galaxies | 3,827 rows, 225 galaxies in `SPARC-rotational.txt` | different file |

For comparison, MOND with one fixed $a_0=1.2\times10^{-10}$ m/s² on the same inputs gives $R^2=0.93$, RMS 23 km/s
and a lower median per-galaxy RMS (10.8 vs 16.7 km/s). These statements concern this file and these choices (mass from
`V_bar(R_last)`, the engine's own surface-density estimate, no external field), not the theory; a different mass
estimate could change them. See `validation/sparc_check.png`.

**CMB peaks.** Planck binned TT gives peaks at $\ell=221,\ 535,\ 814,\ 1123,\ 1409$
(`validation/check_cmb_peaks.py`). The law $\ell_n=279.6\,n$ misses the first peak by 26% and the second and third by
4.5% and 3%, so "within 1%, $A$ fitted to the first peak" fails. A two-parameter line $296.3n-68.4$ is within 3%.

---

## Claims and status

| Claim | Status | Test |
|---|---|---|
| $\kappa_j$ formula, $\kappa_7=1233$, engine constants | verified | `tests/test_claims.py` |
| $\lvert\alpha_G\rvert^2=c\hbar/2GMm$, $=1/2$ at the Planck mass | verified | `tests/test_claims.py` |
| Schwarzschild, Reissner-Nordström, Nariai from the recipe | verified | `tests/test_claims.py`, `tests/test_gaps.py` |
| Exact N-body metric $\sum_a\sigma_a$ | **fails** (not vacuum) | `tests/test_gaps.py` |
| Quantum term resolves the singularity with positive energy | **fails** | `tests/test_gaps.py` |
| Merger at $d=4r_s$ | recorded (depends on which condition is called the merger; total $\sigma=1$ gives $8r_s$) | `tests/test_gaps.py` |
| Constant entropy correction $2\kappa\ell_Q^2/G\hbar$ | **fails** as stated (below) | `tests/test_nariai_entropy.py` |
| SPARC RMS 8.3 km/s, $\chi^2_\nu\approx1.2$, zero free parameters | **fails** | `tests/test_sparc.py` |
| CMB peaks within 1% | **fails** | `tests/test_gaps.py` |
| $\ell_Q=\sqrt{G\hbar^2/c^4}$ is a length; $\rho_{\rm grav}=GM/4c^2r^3$ is an energy density | **fails** in SI | `tests/test_dimensions.py` |
| Neutron-star mass shift about $0.01M_\odot$; fourth-order equation reduces to Einstein; Kerr, Vaidya, interior from the recipe | untested | — |

Predictions listed by the author (none yet confirmed or excluded here): quantum-TOV neutron-star shift (NICER),
equal-mass merger separation, large-scale $\kappa_5$ activation (DESI, Euclid), CMB higher-peak modulation (CMB-S4),
maximum acceleration $a_{\max}=3mc^3/\hbar$; the quantum perihelion ($\sim10^{-90}$ arcsec/century) and GW phase
($\sim10^{-73}$ rad) shifts are unmeasurable. The neutron-star prediction is the main near-term test, and it is untested here.

---

## What the tests found

- **Singularity.** Each term of $f$ has an implied energy density. Mass is vacuum, charge positive, and the quantum
  term $+\kappa\hbar^2/(M^2c^2r^2)$ is the Reissner-Nordström term with $Q^2<0$: **negative** density. The Kretschmann
  scalar still diverges like $r^{-8}$. A regular core needs $f\sim r^2$ at small $r$ (a Hayward-type $f$). Whether the
  fourth-order $\ell_Q^2$ term supplies that is not decided by the static $f$.
- **Entropy** (`validation/nariai_entropy.py`, Planck units). The Euclidean action gives $S_{dS}=3\pi/\Lambda$
  ($S^4$, Euler characteristic 2) and $S_{\rm Nariai}=2\pi/\Lambda$ ($S^2\times S^2$, characteristic 4), matching the
  horizon areas. A truly topological constant is proportional to the Euler characteristic, i.e. to the number of
  horizons: it drops out of $dS/dM$ but not out of de Sitter→Nariai pair creation, where it multiplies the rate by
  $e^{C}$. The documented constant is not that, in either reading of the metric term $A/r^2$:
  - $A=\kappa\ell_Q^2$ constant: the area law gives $4\pi M^2+2\pi A$ (a constant, but $\pi$ times the documented one),
    puts a horizon at $r=\sqrt A$ even at $M=0$, and shifts the black-hole horizon by $+\pi A$ and the cosmological
    one by $-\pi A$, so it is not the same on every horizon.
  - $A=\kappa\hbar^2/(M^2c^2)$ as written in $f$: the shift is $2\pi\kappa/M^2$ and $dM=T\,dS$ fails.
  - A Gauss-Bonnet coupling $\alpha\,\mathrm{GB}$ gives $S=2\pi\alpha\chi$ (checked two ways) but is topological in four
    dimensions and cannot generate the $1/r^2$ term, so the entropy constant and the metric term are independent inputs.
    The entropy table and the metric term need one consistent definition of $A$.
- **Nariai with the quantum term.** The double-root point moves to $r_N=1/\sqrt\Lambda+A/2$, $GM_N=1/(3\sqrt\Lambda)-A/2$.

**Where the checks came from.** The CICY / pyCICY project keeps an exactness ledger and re-derives gravity-side
results by independent routes, including the Nariai horizon $dS_2\times S^2$. Used here: the ledger discipline, the
independent-route principle (the geometry is redone in sympy, not through the QGD code), and the Nariai geometry as the
check on the de Sitter recipe. Nothing in this repository imports CICY.

Run everything with `python -m pytest -q`; expected failures are `xfail` with the reason, so a clean run does not hide
them. Regenerate the ledger with `python validation/claims_ledger.py`.

---

## Open decisions for the author

Fixing these needs a decision about the theory or files not in this repository.

1. **Reconstructed files.** `MasterMetric` and `two_and_three_body_solutions` were missing and are rebuilt from how
   `core/QGD_vs_GR.py` uses them. The script passes 32 of 33 checks; the failure is its own `f_peak` line (GR ISCO
   frequency 67 Hz against a merger peak near 150 Hz). Notably, "merger at $d=4r_s$" is taken as the separation where the
   self part of $\sigma^2$ at the centre of mass equals 1. If you have the originals, they should replace these.
2. **Signature.** `master_metric_diagonal` and `improved_master_metric` give four negative diagonal entries while the
   header declares $(+,-,-,-)$ and $g_{tt}=-(1-f)$. The "improved" $g_{rr}=-(1-f/2)^2(1+f/2)^2$ mixes coordinates, so it
   does not equal isotropic Schwarzschild (`verify_classical_limit` returns `match = False`).
3. **Units.** In SI, $\ell_Q=\sqrt{G\hbar^2/c^4}$ has dimension (mass·length³)$^{1/2}$ (the Planck length is
   $\sqrt{G\hbar/c^3}$), and $\rho_{\rm grav}=GM/4c^2r^3$ has dimension 1/length². Both may be natural-unit shorthand.
4. **Definition of $A$** in the entropy correction (above), and the factor $\pi$ between $2\pi A$ and $2\kappa\ell_Q^2$.
5. **Duplicates and placeholders.** `docs/Ch9-QG_QFT.tex` equals `docs/Ch8-GravitonApplications.tex` (and the PDFs);
   `notebooks/` holds four identical files; `DarkMatter/FullStressEnergyTensor.py` is LaTeX; `core/dark_matter.py` equals
   `DarkMatter/bullet_cluster.py`, as do the two `EFE_solutions_from_QGD_perspective.py`; the citation has placeholders.

Documentation fixed in this revision: $\kappa_7$ (1245→1233); Earth-electron $\lvert\alpha_G\rvert^2=1.46\times10^{39}$
→ $4.36\times10^{-11}$ (classical regime) with electron-electron $2.9\times10^{44}$ as the quantum example, in
`key-equations.md`, `GR-vs-QGD.md`, `docs/complete_paper.tex`, `docs/Ch1-Foundations.tex` and `docs/corrections.md`; broken
README links and the Quick Start; `EFE_solutions_from_QGD_perspective.py` crashing on import; a `SyntaxWarning` in
`core/PN.py`.

---

## Repository structure

```
QGD/
├── DarkMatter/
│   ├── QGD.py
│   ├── dark_matter.py
│   ├── bullet_cluster.py            (identical to core/dark_matter.py)
│   ├── FullStressEnergyTensor.py    (LaTeX source with a .py extension, not Python)
│   ├── Uniqueness_of_kappa_values.py
│   ├── kappa_inversion.py
│   ├── darkmatter-theory.tex, kappa_uniqueness.tex
│   ├── insights.md, kappa-inversion-general.md
│   └── data/                        (SPARC-rotational.txt, PLANCK-TT-Spectrum.txt)
├── core/
│   ├── graviton_field.py, master_metric.py, master-metric2.py
│   ├── qgd_cosmology.py, qgd_energy.py, waveform.py
│   ├── PN.py, PN-generator.py, post-newtonian.py, PostNewtonian/
│   ├── QGD_superposition.py, QGD_vs_GR.py, two_and_three_body_solutions.py (reconstructed)
│   ├── Rosseta_Stone.py, QGD_bug_fixes.py, dark_matter.py
│   └── nrpy/InitialData_QGD.py
├── comparison/
│   ├── EFE_solutions_from_QGD_perspective.py   (also at the repository root)
│   ├── comparison.md, comparison.py
├── validation/
│   ├── dark_matter.py               (QGD engine, v1.8)
│   ├── reproduce_sparc.py, sparc_check.png
│   ├── check_cmb_peaks.py, cmb_peaks_check.png
│   ├── nariai_entropy.py, nariai_entropy.png
│   ├── claims_ledger.py             (writes docs/claims_ledger.md)
│   └── cmb_peaks/
├── docs/                            (theory, LaTeX + PDF)
│   ├── Ch1-Foundations ... Ch16-Quantization (.tex, most with .pdf)
│   ├── complete_paper.tex, corrections.md, claims_ledger.md, summary.md, extended_summary.md, references.bib
├── notebooks/                       (four files with identical content)
├── ongoing-work/                    (drafts: DarkEnergy, inflation, field equations, EOB-correction.py ...)
├── tests/                           (test_claims.py, test_gaps.py, test_dimensions.py, test_sparc.py, test_nariai_entropy.py, test_repo.py, binary_black_hole.py)
├── key-equations.md, GR-vs-QGD.md, comparison-to-GR-and-NG.md, insights-of-QGD.md, QGD_showcase.md
└── requirements.txt
```

---

## Quick Start

```bash
git clone https://github.com/brentharts/Quantum-Gravity-Dynamics.git
cd Quantum-Gravity-Dynamics
pip install -r requirements.txt pandas pytest
python -m pytest -q
```

### Reconstruct the GR solutions from σ-fields

```bash
python comparison/EFE_solutions_from_QGD_perspective.py      # prints "Passed 10/10 checks"
```

### Compute a galaxy rotation curve

```python
import sys
import numpy as np
sys.path.insert(0, "validation")
from dark_matter import predict_rotation_curve

r = np.linspace(1.0, 30.0, 30)                       # kpc
v_gas, v_disk, v_bulge = 20 + 0 * r, 120 * np.tanh(r / 4), 0 * r   # km/s (toy inputs)
v_qgd, kappa = predict_rotation_curve(r, v_gas, v_disk, v_bulge, total_mass=5e10)
print(np.round(v_qgd[::6], 1), np.round(kappa[::6], 2))
```

### Run the waveform module

```bash
python core/waveform.py
```

### Re-run the independent checks

```bash
python validation/reproduce_sparc.py      # SPARC comparison, writes validation/sparc_check.png
python validation/check_cmb_peaks.py      # CMB peak positions
python validation/nariai_entropy.py       # Euclidean vs area entropy, quantum-term shifts
python validation/claims_ledger.py        # regenerate docs/claims_ledger.md
```

## Key References

- Full theoretical derivation: [`docs/complete_paper.tex`](docs/complete_paper.tex)
- Original `.tex` manuscript: available on request
- SPARC rotation curve database: [SPARC](http://astroweb.cwru.edu/SPARC/)
- LIGO GW150914: Abbott et al. (2016), PRL 116, 061102

---

## Citation

```bibtex
@misc{QGD2025,
  title  = {Quantum Gravitational Dynamics: Emergent Geometry from the Dirac Equation},
  author = {[Romeo Matshabba]},
  year   = {2026},
  note   = {GitHub: https://github.com/[matshaba]/Quantum-Gravity-Dynamics}
}
```
