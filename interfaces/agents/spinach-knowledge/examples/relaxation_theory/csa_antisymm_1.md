# examples/relaxation_theory/csa_antisymm_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/csa_antisymm_1.m`
- Signature: `csa_antisymm_1()`
- Total lines: 49

## Purpose

Longitudinal and transverse relaxation rates in a system with a significant antisymmetry in the shielding tensor. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: System specification; implemented by `sys.magnet=14.1`.
- Lines 17-18: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 23-24: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Spinach relaxation rates; implemented by `R=relaxation(spin_system)`.
- Lines 38-40: Textbook relaxation rates; implemented by `[R1Book,R2Book]=rlx_csa(sys.magnet,sys.isotopes{1}, inter.zeeman.matrix{1},inter.tau_c{1})`.
- Lines 42-43: Summary; implemented by `disp(['R1 rate, Spinach: ' num2str(R1Sp)])`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'13C'}`.
- Lines 13: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={[100 20 15`.
- Lines 18: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 19: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 20: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 21: computes `inter.tau_c` using `inter.tau_c={50e-12}`.
- Lines 24: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `R` using `R=relaxation(spin_system)`.
- Lines 33: computes `Lz` using `Lz=state(spin_system,{'Lz'},{1})`.
- Lines 34: computes `Lp` using `Lp=state(spin_system,{'L+'},{1})`.
- Lines 35: computes `R1Sp` using `R1Sp=-(Lz'*R*Lz)/(Lz'*Lz)`.
- Lines 36: computes `R2Sp` using `R2Sp=-(Lp'*R*Lp)/(Lp'*Lp)`.
- Lines 39-40: computes `[R1Book,R2Book]` using `[R1Book,R2Book]=rlx_csa(sys.magnet,sys.isotopes{1}, inter.zeeman.matrix{1},inter.tau_c{1})`.

## Implementation structure

- Longitudinal and transverse relaxation rates in a system
- with a significant antisymmetry in the shielding tensor.
- Calculation time: seconds
- System specification
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Spinach relaxation rates
- Textbook relaxation rates
- Summary

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `state()`, `rlx_csa()`, `num2str()`.
