# examples/relaxation_theory/trosy_selenium.m

- Source: [examples/relaxation_theory/trosy_selenium.m](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/trosy_selenium.m)
- Signature: `trosy_selenium()`

## Purpose

Calculate transverse relaxation matrix elements versus field for selenium and its directly bonded carbon in ethylselenol. The source estimates a runtime of minutes.

## Model and relaxation pathway

The two spins are `77Se` and `13C`. The code reads `../standard_systems/ethylselenol.out`, converts its carbon and selenium data with `g2spinach`, and takes Zeeman matrices and coordinates from converted-data entries 3 (selenium) and 2 (carbon). The source gives conversion parameters [186.38, 0.0] but does not state units for those parameters, tensor values, or coordinates. Both selected nuclei are spin-half isotopes, and this model contains no quadrupolar nucleus or quadrupolar interaction.

Relaxation is explicitly Redfield, with lab-frame relaxation, zero equilibrium, and a 25 ns correlation time; the basis uses `sphten-liouv` with no approximation. The coordinates and shielding tensors feed the dipolar and anisotropic-shielding relaxation terms. The script evaluates normalised selenium and carbon transverse coherences plus opposite-sign two-spin operator combinations (selenium raising coherence paired with selenium-plus-carbon-longitudinal coherence, and the analogous carbon pair). These branches give the TROSY-style relaxation-interference comparison. The source does not separately output a CSA–dipolar cross-correlation term and does not call a stochastic-Liouville solver.

## Calculation and output

Twenty proton Larmor frequencies from 200 to 800 MHz are converted to fields with `2*pi*lin_freq*1e6/spin('1H')`. At every field the spin system, basis, and Redfield relaxation superoperator are rebuilt, and six operator matrix elements are evaluated. Separate plots show selenium and carbon matrix-element branches against proton Larmor frequency in MHz; the vertical axis is relaxation matrix element in Hz. This is a model calculation, not a pulse sequence or measured spectrum.
