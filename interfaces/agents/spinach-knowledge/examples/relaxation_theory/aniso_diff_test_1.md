# examples/relaxation_theory/aniso_diff_test_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/aniso_diff_test_1.m`
- Signature: `aniso_diff_test_1()`
- Total lines: 53

## Purpose

Relaxation superoperator calculation for a dipole-coupled two-spin system with an anisotropic rotational diffusion tensor. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnet field (Tesla); implemented by `sys.magnet=2*pi*950.33e6/spin('1H')`.
- Lines 13-14: Isotopes; implemented by `sys.isotopes={'1H','13C'}`.
- Lines 16-17: Chemical shift tensors (ppm); implemented by `inter.zeeman.eigs={[0.0 0.0 0.0]`.
- Lines 22-23: Cartesian coordinates (Angstrom); implemented by `R=euler2dcm(1,2,3)`.
- Lines 27-28: Scalar couplings (Hz); implemented by `inter.coupling.scalar=cell(2,2)`.
- Lines 31-32: Difusion tensor eigenvalues; implemented by `D=[2.16e8 2.35e8 7.45e8]`.
- Lines 34-35: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 40-41: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 44-45: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 48-49: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=2*pi*950.33e6/spin('1H')`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','13C'}`.
- Lines 17: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[0.0 0.0 0.0]`.
- Lines 19: computes `inter.zeeman.euler` using `inter.zeeman.euler={[0.0 0.0 0.0]`.
- Lines 23: computes `R` using `R=euler2dcm(1,2,3)`.
- Lines 24: computes `inter.coordinates` using `inter.coordinates={[1.13 0.00 0.00]*R`.
- Lines 28: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2,2)`.
- Lines 29: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=145.0`.
- Lines 32: computes `D` using `D=[2.16e8 2.35e8 7.45e8]`.
- Lines 35: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 36: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 37: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 38: computes `inter.tau_c` using `inter.tau_c={1./(6*D)}`.
- Lines 41: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 42: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 45: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- Relaxation superoperator calculation for a dipole-coupled two-spin
- system with an anisotropic rotational diffusion tensor.
- Calculation time: seconds
- Magnet field (Tesla)
- Isotopes
- Chemical shift tensors (ppm)
- Cartesian coordinates (Angstrom)
- Scalar couplings (Hz)
- Difusion tensor eigenvalues
- Relaxation theory
- Basis set
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `spin()`, `euler2dcm()`, `create()`, `basis()`, `relaxation()`.
