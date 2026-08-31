# examples/relaxation_theory/csa_csa_xcorr_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/csa_csa_xcorr_1.m`
- Signature: `csa_csa_xcorr_1()`
- Total lines: 38

## Purpose

Complete Bloch-Redfield-Wangsness relaxation superoperator in a system with two anisotropically shielded nuclei. Spinach relaxation theory mo- dule automatically accounts for all cross-correlations (CSA-CSA cross- correlation is present in this case). Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: System specification; implemented by `sys.magnet=14.1`.
- Lines 20-21: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Relaxation superoperator; implemented by `disp(full(relaxation(spin_system)))`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','13C'}`.
- Lines 15: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[7 15 -22]`.
- Lines 17: computes `inter.zeeman.euler` using `inter.zeeman.euler={[pi/5 pi/3 pi/11]`.
- Lines 21: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 22: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 23: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 24: computes `inter.tau_c` using `inter.tau_c={2e-9}`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- Complete Bloch-Redfield-Wangsness relaxation superoperator in a system
- with two anisotropically shielded nuclei. Spinach relaxation theory mo-
- dule automatically accounts for all cross-correlations (CSA-CSA cross-
- correlation is present in this case).
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Relaxation superoperator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`.
