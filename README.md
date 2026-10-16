# Quantum Gravitational Dynamics (QGD)

> **Core thesis:** Spacetime geometry is not fundamental. It emerges algebraically from a dimensionless phase field $\sigma_\mu$ derived from the semiclassical WKB limit of the Dirac spinor. All of classical General Relativity — and its quantum corrections — follow as consequences.

[![Theory](https://img.shields.io/badge/Theory-QGD-blue)](docs/complete_paper.tex)
[![Solutions](https://img.shields.io/badge/GR_Solutions-Algebraic-green)](comparison/EFE_solutions_from_QGD_perspective.py)
[![Dark Matter](https://img.shields.io/badge/Dark_Matter-No_Per--Galaxy_Fit-orange)](validation/)
[![License](https://img.shields.io/badge/License-MIT-lightgrey)](LICENSE)

---

## Table of Contents

1. [The Core Idea](#the-core-idea)
2. [The Fundamental Chain](#the-fundamental-chain)
3. [Master Equations](#master-equations)
4. [What QGD Solves](#what-qgd-solves)
5. [All GR Solutions Recovered Algebraically](#all-gr-solutions-recovered-algebraically)
6. [Dark Matter as Quantum Structure](#dark-matter-as-quantum-structure)
7. [Falsifiable Predictions](#falsifiable-predictions)
8. [Repository Structure](#repository-structure)
9. [Quick Start](#quick-start)

---

## The Core Idea

General Relativity takes the metric $g_{\mu\nu}$ as fundamental and derives dynamics from the Einstein-Hilbert action. QGD makes a different choice:

**The metric is not fundamental. It is an algebraic output of an underlying phase field.**

Starting from the Dirac spinor $\psi = R(x) e^{iS(x)/\hbar}$ in the semiclassical (WKB) limit, we identify the dimensionless phase gradient:

$$\sigma_\mu(x) \equiv \frac{1}{c} \partial_\mu S(x)$$

This single four-vector encodes all gravitational physics. The metric is then **constructed**, not solved for:

$$\boxed{g_{\mu\nu}(x) = T^\alpha_\mu T^\beta_\nu \left( M_{\alpha\beta} \circ \left[ \eta_{\alpha\beta} - \sum_{a=1}^{N} \varepsilon_a \, \sigma_\alpha^{(a)} \sigma_\beta^{(a)} - \kappa \ell_Q^2 \, \partial_\alpha \sigma^\gamma \partial_\beta \sigma_\gamma \right] \right)}$$

where:
- $T^\alpha_\mu$: coordinate transformation matrix (e.g., spherical, cosmological)
- $M_{\alpha\beta}$: geometric scaling matrix
- $\varepsilon_a \in \{+1, -1\}$: source signature (attractive vs repulsive)
- $\ell_Q = \sqrt{G\hbar^2/c^4}$: quantum gravitational length scale
- $\kappa \approx 2$: quantum stiffness coefficient

---

## The Fundamental Chain

```
Dirac spinor ψ = R(x)·exp(iS/ℏ)
        ↓  WKB limit
Phase field  σ_μ = (1/c)∂_μS
        ↓  algebraic construction
Metric  g_μν = η_μν - Σ εₐ σ_μ^(a) σ_ν^(a) - κℓ_Q² ∂σ∂σ
        ↓  variational principle
Field equation  □_g σ_μ = Q_μ + G_μ + T_μ + κℓ_Q²□_g²σ_μ
        ↓  classical limit ℓ_Q → 0
Einstein's equations  G_μν = (8πG/c⁴) T_μν
```

The relationship to GR is an **equivalence under variable substitution**, not a derivation. When $\sigma_\mu$ satisfies the QGD field equation and $g_{\mu\nu}$ is built from the master metric formula, Einstein's equations are automatically satisfied.

---

## Master Equations

### The Action

$$S_\sigma = \int d^4x \sqrt{-g(\sigma)} \left[ -\frac{c^4}{16\pi G} R[g(\sigma)] + \frac{1}{2}\nabla_\mu \sigma_\nu \nabla^\mu \sigma^\nu - \frac{\ell_Q^2}{2} \nabla_\alpha\nabla_\beta\sigma_\mu \nabla^\alpha\nabla^\beta\sigma^\mu \right] + S_{\text{matter}}$$

### The Field Equation

Variation of $S_\sigma$ with respect to $\sigma_\mu$ yields the fourth-order equation:

$$\boxed{\Box_g \sigma_\mu = Q_\mu(\sigma, \partial\sigma) + G_\mu(\sigma, \ell, H, q) + T_\mu + \kappa\ell_Q^2 \Box_g^2 \sigma_\mu + \mathcal{O}(\ell_Q^4)}$$

Sources:
- $Q_\mu$: nonlinear self-interactions (origin of dark matter phenomenology)
- $G_\mu$: coupling to Kerr-Schild and radiative sectors
- $T_\mu = \frac{1}{2}T^{\mu\nu}\sigma_\nu$: matter stress-energy
- $\kappa\ell_Q^2 \Box_g^2 \sigma_\mu$: quantum gravitational stiffness (resolves singularities)

### The Complete Wavefunction

The four-component gravitational wavefunction encoding all configurations:

$$\psi = \frac{2GMmi}{c\hbar} \left[ \psi_0 \begin{pmatrix}1\\0\\\sqrt{f}\\ ia\sin\theta\sqrt{g}\end{pmatrix} e^{-iS/\hbar} + \psi_1 \begin{pmatrix}0\\1\\-ia\sin\theta\sqrt{g}\\-\sqrt{f}\end{pmatrix} e^{-iS/\hbar} + \psi_2 \begin{pmatrix}\sqrt{f}\\ia\sin\theta\sqrt{g}\\1\\0\end{pmatrix} e^{+iS/\hbar} + \psi_3 \begin{pmatrix}-ia\sin\theta\sqrt{g}\\-\sqrt{f}\\0\\1\end{pmatrix} e^{+iS/\hbar} \right]$$

where the **universal gravitational scalar** $f(r,\theta)$ encodes all physics:

$$\boxed{f(r,\theta) = \underbrace{\frac{2GM}{c^2 r}}_{\text{mass}} - \underbrace{\frac{GQ^2}{c^4 r^2}}_{\text{charge}} + \underbrace{\frac{2Mr}{\Sigma}}_{\text{spin}} + \underbrace{\frac{\Lambda r^2}{3}}_{\Lambda} + \underbrace{\frac{b(r)}{r}}_{\text{pressure}} + \underbrace{H^2(t)r^2}_{\text{expansion}} - \underbrace{\int\frac{P(r)}{\rho(r)c^2}dr}_{\text{EOS}} + \underbrace{\kappa\frac{\hbar^2}{M^2c^2r^2}}_{\text{quantum}}}$$

Every known spacetime is a special case of $f(r,\theta)$.

---

## What QGD Solves

| Long-standing problem | GR status | QGD resolution |
|---|---|---|
| Gravitational energy localization | Pseudotensor only (109-year open problem) | True tensor $T^{\mu\nu}_\sigma$, positive definite |
| Black hole singularities | Generic, unavoidable | Resolved at $r \sim \lambda_C = \hbar/Mc$ |
| N-body exact solutions | No closed form | Exact algebraic: $\sigma_t = \sum_a \sqrt{2GM_a/c^2 r_a}$ |
| Dark matter | Requires new particles | Factorial $\kappa_j$ structure of quantum phase Taylor expansion |
| Quantum corrections | Undefined in GR | Explicit $\mathcal{O}(\hbar^2)$: $\delta g_{tt} = -G\hbar^2/(Mc^4 r^3)$ |
| Binary waveforms | Supercomputer, weeks | Algebraic, O(N²), seconds |
| Cosmological constant | Fine-tuning problem | $\rho_\sigma = 3H_0^2/8\pi G$, $w = -1$ from attractor solution |

---

## All GR Solutions Recovered Algebraically

The **source signature recipe** — choose $\sigma^{(a)}$ and $\varepsilon_a$, apply the master metric:

| Spacetime | $\sigma$-field | $\varepsilon_a$ | Physical effect |
|---|---|---|---|
| Schwarzschild | $\sqrt{2GM/c^2r}$ | $+1$ | Attractive mass |
| Kerr | $a\sin\theta\sqrt{2GM/c^2r}$ | $+1$ | Frame dragging |
| Reissner-Nordström | $\sqrt{GQ^2/c^4r^2}$ | $-1$ | EM repulsion |
| de Sitter ($\Lambda > 0$) | $Hr$ | $+1$ | Cosmological expansion |
| Anti-de Sitter | $\|H\|r$ | $-1$ | AdS geometry |

**Frame dragging as interference:** The Kerr off-diagonal term $g_{t\phi}$ is not input — it emerges automatically from the cross-product $\sigma_t^{(\text{mass})} \times \sigma_\phi^{(\text{spin})}$.

See [`comparison/EFE_solutions_from_QGD_perspective.py`](comparison/EFE_solutions_from_QGD_perspective.py) for the algebraic constructions of each metric (it checks 10 of them and prints the result).

---

## Dark Matter as Quantum Structure

The Taylor expansion of the gravitational phase factor:

$$e^{i\phi} = e^{i\sigma_\mu x^\mu/\hbar} = \sum_{j=0}^{\infty} \frac{(i\sigma)^j}{j!}$$

generates factorial enhancement factors:

$$\kappa_j = \sqrt{\frac{(2j-1)!}{2^{2j-2}}} \qquad \Rightarrow \qquad \kappa = [1.00,\ 1.225,\ 2.74,\ 8.87,\ 37.7,\ 197,\ 1233, \ldots]$$

**Modified rotation curve (no per-galaxy fit parameters; the engine's shared empirical constants are listed under "Independent checks"):**

$$v^2(r) = \frac{GM_{\text{baryon}}(<r)}{r}\left[1 + \sum_{j=2}^{4} \kappa_j \, g_j(r)\right]$$

**Validation reported by the author (4,248 measurements, 175 galaxies; not reproduced in full, see "Independent checks"):**
- $R^2 = 0.908$ (reproduced on the shipped data file: 0.91)
- $\chi^2_\nu \approx 1.2$ (not reproduced: about 150 under the choices in `validation/reproduce_sparc.py`)
- RMS = 8.3 km/s (not reproduced: about 26 km/s)
- No per-galaxy fit parameters (vs 5–7 for $\Lambda$CDM); the shared empirical constants of the engine are not counted in that comparison

Cross-dataset universality: same $\kappa_j$ values for clusters, ellipticals, CMB peak spacing, and wide binary External Field Effect.

See [`validation/`](validation/) and [`validation/dark_matter.py`](validation/dark_matter.py).

---

## Falsifiable Predictions

All predictions follow from the theory with no additional assumptions:

| Prediction | Value | Testable with |
|---|---|---|
| Neutron star mass shift | $\sim 0.01 M_\odot$ from quantum TOV | NICER (current) |
| Binary merger separation | $d = 4r_s$ for equal masses | Next-gen GW detectors |
| Large-scale $\kappa_5$ activation | Specific correlation at 10–100 Mpc | DESI, Euclid |
| CMB higher peak modulation | $\ell_n \propto \kappa_{j(n)} \cdot f(n)$ | CMB-S4 |
| Maximum acceleration | $a_{\max} = 3mc^3/\hbar$ | — |
| Quantum perihelion shift | $\sim 10^{-90}$ arcsec/century | Unmeasurable |
| GW phase quantum shift | $\sim 10^{-73}$ rad | Unmeasurable |

The neutron star mass prediction with NICER is the **critical near-term falsification test**.

---

## Independent checks

`tests/` and `validation/reproduce_sparc.py` re-check some of the statements above without importing the
metric code (the algebra is redone in sympy).  Run `pip install -r requirements.txt pytest pandas`, then
`python -m pytest -q`.  Expected failures are marked `xfail` with the reason, so a clean run does not hide them.

**Passing:** the factorial formula for the `kappa_j` and the constants in `validation/dark_matter.py`;
`alpha_G` identities and `|alpha_G|^2 = 1/2` at the Planck mass; the Kerr-Schild form of the Schwarzschild
recipe is Ricci-flat and the Reissner-Nordstrom recipe (`epsilon = -1` for charge) has exactly the Maxwell
stress, while the opposite sign gives the opposite-sign stress.

**Fixed in this revision:**
- `kappa_7` in the README (1245 -> 1233).
- The Earth-electron value `|alpha_G|^2 = 1.46e39` in `key-equations.md`, `GR-vs-QGD.md` and `docs/complete_paper.tex`
  is now 4.36e-11 (classical regime), and electron-electron 2.9e44 is given as the quantum-regime example;
  `docs/Ch1-Foundations.tex` had 4.35e-10 (a factor 10 off) and now reads 4.36e-11; `docs/corrections.md` no longer
  says the electron-electron value is ~1e-45 (that number is `G m_e^2/(hbar c)`).
- README links to files that do not exist, the repository-structure listing, and the Quick Start (it imported
  modules that are not in the repository); the Quick Start is now executed by `tests/test_repo.py`.
- `EFE_solutions_from_QGD_perspective.py` (both copies) crashed on import; it now runs and passes its 10 checks.
- A `SyntaxWarning` in `core/PN.py`, and the `master_metric_diagonal` docstring (it returns only the matrix).

**Not reproduced:**
- `SPARC-rotational.txt` holds 3,827 rows and 225 galaxies; the README quotes 4,248 measurements and 175 galaxies.
- On that file, with the choices listed in `reproduce_sparc.py` (mass from `V_bar(R_last)`, engine's own surface-density
  estimate, no external field) the engine gives pooled R^2 = 0.91, matching the README, but RMS = 26 km/s and
  chi^2_nu = 152, not 8.3 km/s and 1.2.  MOND with one fixed `a0 = 1.2e-10 m/s^2` and the same inputs gives
  R^2 = 0.93, RMS = 23 km/s and a lower median per-galaxy RMS (10.8 vs 16.7 km/s).
- The engine has empirical parameters (`SIGMA_CRIT`, `ALPHA`, `BETA_0`, `LOG_M_TRIGGER`, the two mass-to-light
  ratios, `G_CRIT`), so "zero free parameters per galaxy" holds only in the sense that they are shared across galaxies.
  Refitting the first three on a random half of the galaxies improves the held-out half only slightly
  (RMS 24.8 vs 25.6 km/s in the first split), and the fit runs to its bounds, so they are poorly constrained here.

These are statements about this file and these choices, not about the theory; a different mass estimate or
subset could change the numbers.  `validation/sparc_check.png` shows the comparison.

## Known issues

Found by the checks above and left for the author, because fixing them needs a decision about the theory or
files that are not in this repository:

1. `core/QGD_vs_GR.py` needed `MasterMetric` and `two_and_three_body_solutions`, which were missing.  Both are now
   **reconstructed** from how the script uses them (`core/master_metric.py`, `core/two_and_three_body_solutions.py`);
   the script runs and passes 32 of 33 checks.  The failing one is its own `f_peak` comparison: it computes the GR ISCO
   frequency of the total mass (67 Hz) and compares it with the ~150 Hz merger peak.  Two reconstruction choices are
   documented in the code, notably that the "merger at d = 4 r_s" is the separation where the *self* part of sigma^2
   at the centre of mass equals 1 (the *total* sigma = 1 would give 8 r_s).  If you have the original files, they
   should replace these.
2. `core/master_metric.py`: `master_metric_diagonal` and `improved_master_metric` give a metric with all four
   diagonal entries negative, i.e. no Lorentzian signature, while the header declares (+,-,-,-) and `g_tt = -(1-f)`.
   The "improved" `g_rr = -(1-f/2)^2 (1+f/2)^2` mixes coordinates: isotropic Schwarzschild has
   `g_tt = -(1-f/2)^2/(1+f/2)^2` and `g_rho,rho = (1+f/2)^4` with `f = r_s/(2 rho)`, so the claim that it matches
   isotropic Schwarzschild exactly does not hold (`verify_classical_limit` returns `match = False`).
3. `docs/Ch9-QG_QFT.tex` is byte-for-byte `docs/Ch8-GravitonApplications.tex` (and the two PDFs are identical).
4. `notebooks/` holds four files with identical content (one has no `.ipynb` extension).
5. `DarkMatter/FullStressEnergyTensor.py` is a LaTeX document with a `.py` extension.
6. `core/dark_matter.py` and `DarkMatter/bullet_cluster.py` are identical, as are the two copies of
   `EFE_solutions_from_QGD_perspective.py`.
7. Unit problems in the formulas as written (`tests/test_dimensions.py`): `ell_Q = sqrt(G hbar^2 / c^4)` has dimension
   (mass x length^3)^(1/2), not length (the Planck length is `sqrt(G hbar / c^3)`); `rho_grav = G M / (4 c^2 r^3)` has
   dimension 1/length^2, not energy density.  Both may be natural-unit shorthand; the SI forms need fixing.
8. The citation block still has placeholders (`[author]`, `[matshaba]`) and the BibTeX key says 2025 while the year is 2026.

---

## Where the theory has gaps (tests and ledger)

`tests/test_gaps.py` turns suspected gaps into tests and `docs/claims_ledger.md` (generated by
`validation/claims_ledger.py`) lists every claim checked, with its status.  What the tests show:

- **Superposition.** The recipe `g = eta + sum eps sigma sigma` is exact when the terms share one null direction
  (mass + charge on one centre: the mixed Ricci tensor of a Kerr-Schild metric is linear in `h`), which is why the
  Schwarzschild, Reissner-Nordstrom and de Sitter rows work.  It is **not** exact for two separated masses: the
  `Sigma_a sigma_a` N-body metric has a nonzero Ricci tensor, of order `(GM)^2` (the cross term).  So the "exact N-body
  solution" is a statement about the sigma-field sum, not about a vacuum solution of Einstein's equations.
- **Singularity resolution.** In the static form the universal `f` is additive in the Einstein tensor, so each term
  has a definite implied energy density: the mass term is vacuum, the charge term positive, and the quantum term
  `+kappa hbar^2/(M^2 c^2 r^2)` is the Reissner-Nordstrom term with `Q^2 < 0`, i.e. **negative** density, and the
  Kretschmann scalar still diverges like `r^-8`.  A regular core needs `f ~ r^2` at small `r` (e.g. a Hayward-type
  `f`).  Whether the fourth-order `l_Q^2` term supplies that is not decided by the static `f`.
- **de Sitter and Nariai.** The `sigma_Lambda` recipe reproduces the Schwarzschild-de Sitter double-horizon (Nariai)
  limit exactly: `r = 1/sqrt(Lambda)`, area `4 pi / Lambda`, `G M = c^3 / (3 sqrt(3) H)`.
- **Entropy** (`validation/nariai_entropy.py`, `tests/test_nariai_entropy.py`, Planck units).  The Euclidean action
  gives `S_dS = 3 pi/Lambda` (`S^4`, chi = 2) and `S_Nariai = 2 pi/Lambda` (`S^2 x S^2`, chi = 4), matching the horizon
  areas.  A truly topological constant is proportional to chi, i.e. to the number of horizons: it drops out of `dS/dM`
  but **not** out of de Sitter -> Nariai pair creation, where it would multiply the rate by `exp(C)`.  The documented
  constant `2 kappa l_Q^2/(G hbar)` is not that, in either reading of the metric term:
  (I) with `A = kappa l_Q^2` held constant the area law gives `4 pi M^2 + 2 pi A` (a constant, but `pi` times the
  documented one) and a horizon at `r = sqrt(A)` even at `M = 0`; the shift is `+pi A` on black-hole-type horizons and
  `-pi A` on the cosmological one, so it is not the same on every horizon; (II) with `A = kappa hbar^2/(M^2 c^2)` as
  written in the f formula the shift is `2 pi kappa/M^2` and `dM = T dS` fails.  The entropy table and the metric term
  need one consistent definition of `A` before this can be tested further.
- **CMB peaks** (`validation/check_cmb_peaks.py`, Planck binned TT): measured peaks at l = 221, 535, 814, 1123, 1409;
  the law `l_n = 279.6 n` misses the first by 26% and the second and third by 4.5% and 3%, so "within 1%" and
  "A fitted to the first peak" do not hold.  A two-parameter line `l_n = 296.3 n - 68.4` stays within 3%.

**Connection to CICY.** The CICY/pyCICY project keeps an exactness ledger (exact, needs outside input, does not
exist) and re-derives gravity-side results by independent routes, including the Nariai horizon `dS_2 x S^2`.  Used here:
the ledger discipline (`docs/claims_ledger.md`), the independent-route principle (every test above redoes the geometry
with sympy rather than calling the QGD code), and the Nariai geometry as the check on the de Sitter recipe.  Nothing in
the QGD repository imports CICY.

---

## Repository Structure

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
│   ├── reproduce_sparc.py, sparc_check.png, check_cmb_peaks.py, cmb_peaks_check.png, claims_ledger.py
│   └── cmb_peaks/
├── docs/                            (theory, LaTeX + PDF)
│   ├── Ch1-Foundations ... Ch16-Quantization (.tex, most with .pdf)
│   ├── complete_paper.tex, corrections.md, claims_ledger.md, summary.md, extended_summary.md, references.bib
├── notebooks/                       (four files with identical content, see Known issues)
├── ongoing-work/                    (drafts: DarkEnergy, inflation, field equations, EOB-correction.py ...)
├── tests/                           (test_claims.py, test_gaps.py, test_dimensions.py, test_sparc.py, test_repo.py, binary_black_hole.py)
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

The engine has empirical parameters (see "Independent checks"); it is not parameter-free.

### Run the waveform module

```bash
python core/waveform.py
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
