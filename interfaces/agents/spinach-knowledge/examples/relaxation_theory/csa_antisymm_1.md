# examples/relaxation_theory/csa_antisymm_1.m

- Signature: `csa_antisymm_1()`

## Purpose

Calculate longitudinal and transverse relaxation rates for a single `13C` nucleus with a shielding tensor that has a significant antisymmetric component, then compare Spinach projections with textbook CSA rates. Calculation time: seconds.

## Model and parameters

- Field: `14.1 T`. The shielding matrix (ppm) is `[100 20 15; 20 0 30; 25 10 -30]`.
- Redfield relaxation uses `tau_c={50e-12}`, zero equilibrium, and `labframe` retention; the basis is `sphten-liouv` with no approximation.

## Calculation

After computing `R=relaxation(spin_system)`, the example projects `R` onto `Lz` and `L+` to obtain `R1Sp` and `R2Sp`. It calls `rlx_csa` with the same field, isotope, shielding matrix, and correlation time for `R1Book` and `R2Book`, then prints both pairs of values.
