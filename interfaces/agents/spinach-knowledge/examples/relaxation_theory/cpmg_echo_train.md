# examples/relaxation_theory/cpmg_echo_train.m

- Signature: `cpmg_echo_train()`

## Purpose

Simulate and plot a CPMG echo train over a powder grid for a two-proton system. Calculation time: seconds.

## System and relaxation model

- Field: `14.1 T`; isotopes: `{'1H','1H'}`. The shielding principal-value expressions are `[-2 -2 4]-5` and `[-1 -3 4]+5`, with zero Euler angles.
- The `t1_t2` model uses `r1_rates={50.0 50.0}` and `r2_rates={150.0 150.0}`, zero equilibrium, and secular relaxation retention.
- The basis is `sphten-liouv` without approximation; trajectory-level processing is disabled.

## Powder simulation

The powder grid is `rep_2ang_200pts_sph`. The experiment selects `1H`, uses `L+` for the initial state and coil, and sets `Lx` as the pulse operator. It runs `powder(spin_system,@cpmg,parameters,'nmr')` with 10 loops, a `1e-5 s` timestep, and 100 points, then plots the real FID versus time.
