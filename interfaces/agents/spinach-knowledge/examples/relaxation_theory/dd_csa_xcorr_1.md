# examples/relaxation_theory/dd_csa_xcorr_1.m

- Signature: `dd_csa_xcorr_1()`

## Purpose

Construct and display a complete Redfield relaxation superoperator for anisotropically shielded `1H` and `13C` nuclei with a dipolar interaction. The source identifies CSA–CSA and DD–CSA cross-correlations as present and states that the dipolar coupling is computed from the spin coordinates. Calculation time: seconds.

## Model and parameters

- Field: `14.1 T`; shielding principal-value rows (ppm): `[7 15 -22]` and `[11 18 -29]`.
- Shielding Euler-angle rows: `[pi/3 pi/4 pi/5]` and `[pi/6 pi/7 pi/8]`.
- Coordinates (Å): `[0 0 0]` and `[0 0 1.02]`.
- Redfield relaxation uses `tau_c={1e-9}`, zero equilibrium, and `labframe` retention; the basis is `sphten-liouv` with no approximation.

## Calculation

Spinach creates and bases the system, evaluates `relaxation(spin_system)`, and prints the full superoperator.
