# examples/nmr_liquids/roesy_sucrose.m

- Signature: `roesy_sucrose()`

## Purpose

Simulates a ROESY spectrum of sucrose using magnetic parameters computed with DFT. Stated calculation time: minutes.

## Spin system and simulation settings

- Loads proton spin-system data from `../standard_systems/sucrose.log` using `gparse` and `g2spinach`, with `options.min_j=1.0`; sets the magnetic field to `5.9`.
- Uses the `sphten-liouv` formalism, `IK-2` approximation, scalar-coupling connectivity, and proximity level `3`.
- Specifies Redfield relaxation, zero equilibrium, secular relaxation terms, and correlation time `200e-12`.
- Enables `greedy`, **disables `krylov`**, and sets `sys.tols.prox_cutoff=4.0`.
- Creates the spin system and basis, then runs `liquid(spin_system,@roesy,parameters,'nmr')` with mixing time `0.5`, offset `800`, sweep `[1700 1700]`, `512 × 512` points, `2048 × 2048` zero filling, proton spins, ppm axes, and initial state `Lz` for `1H`.

## Processing and output

- Applies squared-cosine apodisation in both dimensions to the cosine and sine FIDs.
- Fourier-transforms the cosine and sine data, combines their imaginary and real components as a States signal, and Fourier-transforms the result to obtain the spectrum.
- Plots `real(spectrum)` with `plot_2d`.