# examples/relaxation_theory/t1t2_strychnine.m

- Signature: `t1t2_strychnine()`

## Purpose

Relaxation analysis for strychnine, dipolar processes only. Calculation time: seconds.

## Physical / mathematical content

- The spin system is obtained from `strychnine({'1H'})`, and the magnetic field is set to `5.9`.
- The relaxation settings are `inter.relaxation={'redfield'}`, `inter.equilibrium='zero'`, `inter.rlx_keep='kite'`, and `inter.tau_c={200e-12}`.

## Numerical / algorithmic content

- The basis uses `sphten-liouv` formalism, `IK-2` approximation, `scalar_couplings` connectivity, and proximity level `3`.
- The distance cut-off is set by `sys.tols.prox_cutoff=4.0`.

## Implementation structure

- The function creates the spin system with `create(sys,inter)`, applies `basis(spin_system,bas)`, and runs `relaxan(spin_system)`.
