# examples/relaxation_theory/t1t2_strychnine.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/examples/relaxation_theory/t1t2_strychnine.m)

## Purpose

This example runs a relaxation analysis for the proton system of strychnine, using dipolar processes only, as stated in the source. It is a calculated model, not an experimental T1/T2 measurement. The source lists a calculation time of seconds.

## Spin system and relaxation model

The system is initialised by `strychnine({'1H'})`, then sets `sys.magnet=5.9`. The basis uses `sphten-liouv`, the `IK-2` approximation, connectivity from scalar couplings, and proximity level 3. Redfield relaxation is selected with zero equilibrium, `kite` relaxation retention, and `inter.tau_c={200e-12}`. A distance cutoff of 4.0 is applied.

## Analysis

The script creates the Spinach system, constructs the basis, and calls `relaxan` for the relaxation analysis. It does not define a pulse sequence, initial state, detection operator, spectral sweep, or custom plot; any reported relaxation quantities come from the analysis routine. No numerical T1/T2 values are hard-coded in the example.
