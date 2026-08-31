# examples/relaxation_theory/t1t2_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/t1t2_strychnine.m`
- Signature: `t1t2_strychnine()`
- Total lines: 39

## Purpose

Relaxation analysis for strychnine, dipolar processes only. Calculation time: seconds.

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Spin system properties; implemented by `[sys,inter]=strychnine({'1H'})`.
- Lines 13-14: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 16-17: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 22-23: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 28-29: Distance cut-off; implemented by `sys.tols.prox_cutoff=4.0`.
- Lines 31-32: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 35-36: Relaxation analysis; implemented by `relaxan(spin_system)`.

### Key state/data transformations

- Lines 11: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H'})`.
- Lines 14: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 17: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 18: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 19: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 20: computes `bas.space_level` using `bas.space_level=3`.
- Lines 23: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 24: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 25: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 26: computes `inter.tau_c` using `inter.tau_c={200e-12}`.
- Lines 29: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 32: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- Relaxation analysis for strychnine, dipolar processes only.
- Calculation time: seconds.
- Spin system properties
- Magnet field
- Basis set
- Relaxation theory parameters
- Distance cut-off
- Spinach housekeeping
- Relaxation analysis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `relaxan()`.
