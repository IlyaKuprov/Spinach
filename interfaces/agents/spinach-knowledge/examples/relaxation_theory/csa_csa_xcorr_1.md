# examples/relaxation_theory/csa_csa_xcorr_1.m

- Signature: `csa_csa_xcorr_1()`

## Purpose

Construct and display the complete Redfield relaxation superoperator for an anisotropically shielded `1H`–`13C` pair. The source notes that the CSA–CSA cross-correlation contribution is present and that the Spinach relaxation module accounts for the cross-correlations. Calculation time: seconds.

## Model and parameters

- Field: `14.1 T`; shielding principal-value rows (ppm): `[7 15 -22]` and `[11 18 -29]`.
- Shielding Euler-angle rows: `[pi/5 pi/3 pi/11]` and `[pi/6 pi/7 pi/15]`.
- Redfield relaxation uses `tau_c={2e-9}`, zero equilibrium, and `labframe` retention. The basis is `sphten-liouv` with no approximation.

## Calculation

The function creates and bases the system, evaluates `relaxation(spin_system)`, and prints the full superoperator. The source defines no dipolar coupling or coordinates for this example.
