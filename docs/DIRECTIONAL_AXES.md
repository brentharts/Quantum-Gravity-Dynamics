# Directional axes and timelines: what fits QGD and what does not

This note brings the directional-axis toy models of the `spectre-monotile-py` work into this repository, where the checking culture already exists: every
statement gets a status and a test (`directional/LEDGER.md`, same statuses as `docs/claims_ledger.md`). Nothing here is a claim that QGD is right or
that the toy models describe nature.

## The models

In `spectre-monotile-py/research/toymodel_v1` a particle history is a list of steps, each one of six signed directions (+X, -X, +Y, -Y, +Z, -Z). A
"consistency contract" edits steps before and after a forced step until the net flux is small; a density profile sets a clock rate `1/rho`; coupled pairs
start anti-correlated; a grid of timelines carries a source. The motivation is a search for combinatorial structures with an arrow of time: "direction"
(how many, which sign) rather than "dimension", and no assumption that space is 3+1.

## Where it fits QGD

- **Time dilation is the one place the two can be compared directly.** QGD gives the clock of a static observer, `sqrt(1 - sigma_t^2)` with
  `sigma_t^2 = 2GM/(c^2 r)` (Kerr-Schild form, redone in sympy in `tests/test_directional.py`). The toy rule `clock = 1/rho` agrees with it to first order if
  `rho = 1 + GM/(c^2 r)` and differs at second order (`directional/clock_rate_vs_sigma.py`). The density that matches exactly is `(1 - r_s/r)^(-1/2)`.
  That is a definition of `rho`, not an emergence of the dilation, and the Lorentzian well of the toy has the wrong tail (`d^-2` against `r^-1`) and no horizon.
- **The checking style.** This repository already keeps the same ledger discipline as the Exact / Numeric / Hypothesis labels of the tile work, so the
  toy models can be put through it without changing the way either project works.

## What does not fit (yet)

- **E8.** QGD is stated in four dimensions with no extra structure; the $E_8$ results (Hopf fibres, the Spectre lift lattice inside $E_8$, the 192 count)
  live in rank 8 and stay in `spectre-monotile-py`. They are cited, not copied, and no step here uses them.
- **Phase field and directions.** `sigma_mu` is a four-component phase gradient, and `tests/test_sigma_phase_gradient.py` records that the recipe is not even
  closed. The six directions are discrete steps. No map between them is given, and the ledger says so (UNTESTED).
- **The coupled network.** The paper draft says it gives entanglement without action at a distance. As written it does the opposite (below).
- **The grid.** `GridWaveSimulator.py` is described as giving travelling waves; not tested yet (UNTESTED).

## What the first tests found

1. **Counts are exact.** Closed timelines (zero net flux) of n steps number 1, 0, 6, 0, 90, 0, 1860, 0, 44730 for n = 0..8. The contract's stopping rule
   `|flux|^2 <= 3` is, by parity, `|flux|^2 in {0, 2}` for an even number of steps.
2. **The solver scores its edits with a sign that looks wrong.** `AxisNetwork._resolve_consistency` minimises `|old + new + flux|^2`; the flux after replacing a
   step is `|flux - old + new|^2`. With the corrected sign the contract reaches the target in 89% of 2000 seeded runs (n = 12) against 66% as written.
   (This is a statement about the code as it stands in `spectre-monotile-py`; the original may have intended something else.)
3. **Time dilation is an input.** See above. If the paper draft says that dilation is *derived* from the interaction density, the supported statement is
   that a prescribed density with `clock = 1/rho` reproduces the QGD / GR clock to first order in `GM/(c^2 r)`.

4. **The contract copies.** In `CoupledParticleNetwork` the node contract returns its reference axis whatever the current axis is, so measuring A along an
   axis overwrites A's past up to the measured step and all of B with the inverse axis. B's history then identifies the axis A measured with certainty
   (log2 3 bits): the model as written signals instantly. The same copy rule in `GravitationalTimelineNetwork` erases the forced event under uniform density.
5. **CHSH cannot separate the toy from quantum mechanics inside the six directions.** There the toy's correlation is -a.b, exactly the singlet's, and the
   maximum CHSH value is 2 for both. If B is read along any direction by the simplest rule (sign of the projection of B's vector), the correlation is -sgn(a.b) and
   CHSH reaches 4 against the singlet's 2 sqrt 2 for the same settings, with each outcome still unbiased. The readout rule is an assumption of ours; the
   toy has none.

6. **Four contract policies, compared exactly** (`directional/contract_options.py`). Copy (as written): B's history is overwritten, CHSH 4, signals completely. Global
   conservation (B repairs the joint flux): signals through B's flux (total variation 0.161), and with a flux readout stays below 2. Local conservation: no signalling and no
   correlation. Revelation of a pre-existing origin axis: local, exact maximum 4/3. So within these policies a model signals, or does not correlate, or is a local hidden-variable
   model. This is not a no-go for retrocausal models in general: a model whose history ensemble depends on both settings (settings as future boundary conditions) is outside
   these policies and is the natural next thing to build.
7. **The derivation that does hold** (`directional/update_rate_clock.py`). If a node updates with probability 1/rho per sweep, its tick rate is 1/rho, and the redshift between two stations is the
   density ratio; with rho = (1 - r_s/r)^(-1/2) it is the Schwarzschild redshift. The density is still an input.
8. **Two sources expose a first-order difference** (`directional/two_source_clocks.py`). The toy with additive density agrees with linearised GR at first order. The QGD sigma-sum
   clock does not: sigma_a is proportional to sqrt(G M_a), so the cross term 2 sigma_1 sigma_2 is of order G sqrt(M_1 M_2), the same order as the Newtonian terms. Far from two masses
   the sigma-sum behaves like the single mass (sqrt M_1 + sqrt M_2)^2 (twice the Newtonian mass for equal masses), and the cross term is not harmonic. Consistently, the Ricci residual of the
   sigma-sum metric halves when the masses halve (first order), against the quartering of the linear Kerr-Schild sum that `tests/test_gaps.py` checks. `README.md` previously said "order (GM)^2" for
   the sigma-sum; that sentence is corrected. Caveat: this is about the formulas as printed in `key-equations.md` section 13; if the sigma_tot of the superposition theorem is meant to be something
   other than the sum, the statement does not apply to it.
9. **Settings as boundary conditions, priced** (`directional/retro_boundary_model.py`). Let the history ensemble depend on both settings, with the origin axis as the hidden
   state, outcomes read from it, and no signalling. A linear program gives the largest CHSH value as `min(4, 4/3 + 8 delta)`, where delta is how far each hidden-state probability may move
   from the uniform 1/6. No dependence (delta = 0) is the local model (4/3); a Bell violation needs delta > 1/12 (a 50% relative change of each probability); 4 is reached at delta = 1/3.
   The natural weight rule "favour the origin axes that make the outcomes anti-correlated" does not pay the price cleanly: it signals for every kappa > 0, and reaches CHSH 2.75 only at signalling 0.89.
   So retro-dependence in this toy can in principle give a Bell violation without signalling, but it costs a large dependence of the ensemble on the settings, and the rule that would
   supply it is not known. (The size of the dependence needed is the quantity the literature on relaxing measurement independence studies; we have not read it, and cite nothing.)

## Next

1. **A rule that pays the price.** The program shows what dependence is needed (delta > 1/12) but no rule. Candidates to test: weights from a path action with boundary terms at both ends and a normalisation over the whole history (a two-state-vector structure), so that the dependence on the two settings is not put in by hand; check no signalling and the CHSH value exactly.
2. **Wave speed on the grid.** Check `GridWaveSimulator.py` against the local clock rate now that the clock is derived from the update rate.
3. **The binary icosahedral group in this repository.** `README.md` describes the Poincare dodecahedral space as S^3 / I*, with `validation/dodecahedral_spectrum.py`; that file
   is not in this clone. The 24-cells inside the 600-cell and their coset partitions found in `spectre-monotile-py` (25 inscribed, 10 partitions; index 5 of 2T in 2I) are
   exactly the structure of the 5-fold cover S^3/2T -> S^3/2I and could be added as an independent check on the multiplicities there.

## Run

    pip install -r requirements.txt pytest
    python directional/axis_timelines.py
    python directional/clock_rate_vs_sigma.py
    python directional/coupled_network.py
    python directional/contract_options.py
    python directional/two_source_clocks.py
    python directional/update_rate_clock.py
    python directional/retro_boundary_model.py        # about 50 s
    python -m pytest -q tests/test_directional.py tests/test_directional_coupled.py tests/test_directional_models.py tests/test_directional_retro.py
