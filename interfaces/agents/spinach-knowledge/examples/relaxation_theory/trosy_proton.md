# examples/relaxation_theory/trosy_proton.m

- Source: [examples/relaxation_theory/trosy_proton.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/trosy_proton.m)
- Signature: `trosy_proton()`

## Purpose

Calculate transverse relaxation matrix elements versus field for the C–H group at position 3 of a tyrosine aromatic ring. The source estimates a runtime of minutes.

## Model and relaxation pathway

This is a two-spin `1H`/`13C` model. The code parses `../standard_systems/amino_acids/tyr.log` as a 3-fluorotyrosine DFT calculation, maps carbon and hydrogen to `13C` and `1H` with conversion arguments `[186.38, 33.44]`, and takes the proton and carbon Zeeman matrices and coordinates from entries 10 and 5 of the converted DFT data, respectively. No unit for the imported tensor values or coordinates is stated in this source. Relaxation is explicitly set to Redfield, with lab-frame relaxation, zero equilibrium, and a 25 ns correlation time; the basis uses `sphten-liouv` with no approximation.

At each field the code forms normalised proton and carbon raising-operator coherences, then their left/right combinations with two-spin terms (proton raising coherence paired with proton-plus-carbon-longitudinal coherence, and carbon raising coherence paired with carbon-plus-proton-longitudinal coherence). The opposite-sign branches compare the relaxation interference used in a TROSY-style analysis. The DFT shielding tensors and the C–H geometry supply the anisotropic-shielding and dipolar interaction information for the Redfield calculation; the source does not report separately resolved CSA–dipolar cross-correlation values. The script calls the Redfield model, not a stochastic-Liouville solver.

## Calculation and output

Twenty proton Larmor frequencies from 200 to 800 MHz are converted to fields using `2*pi*lin_freq*1e6/spin('1H')`. The script rebuilds the spin system and basis and evaluates the relaxation superoperator at each field. It plots the normalised proton and carbon operator matrix elements and their two opposite-sign branches in separate figures. The vertical axis is relaxation matrix element in Hz; the horizontal axis is proton Larmor frequency in MHz. There is no pulse sequence, time-domain acquisition, or measured spectrum in this example.
