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

## Next

1. **A contract that uses the current state.** For example, a node keeps its axis unless a conservation rule forces a change, with B's history revealed
   only through a stated readout rule. Then repeat the CHSH computation; the target is not to exceed 2 sqrt 2 without signalling.
2. **Spatial clocks.** Put the QGD density on a grid and check the grid simulator's wave speed against the local clock.
3. **The binary icosahedral group in this repository.** `README.md` describes the Poincare dodecahedral space as S^3 / I*, with `validation/dodecahedral_spectrum.py`; that file
   is not in this clone. The 24-cells inside the 600-cell and their coset partitions found in `spectre-monotile-py` (25 inscribed, 10 partitions; index 5 of 2T in 2I) are
   exactly the structure of the 5-fold cover S^3/2T -> S^3/2I and could be added as an independent check on the multiplicities there.

## Run

    pip install -r requirements.txt pytest
    python directional/axis_timelines.py
    python directional/clock_rate_vs_sigma.py
    python directional/coupled_network.py
    python -m pytest -q tests/test_directional.py tests/test_directional_coupled.py
