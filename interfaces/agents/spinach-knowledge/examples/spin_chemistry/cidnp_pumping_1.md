# examples/spin_chemistry/cidnp_pumping_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/spin_chemistry/cidnp_pumping_1.m`
- Signature: `cidnp_pumping_1()`
- Total lines: 69

## Purpose

A simulation of the matrix in Equation 2 of IK's paper on chemically amplified NOEs (https://doi.org/10.1016/j.jmr.2004.01.011). Calculation time: seconds.

## Physical / mathematical content

- Spin-chemistry examples. These scripts treat radical pairs, recombination channels, chemically induced dynamic nuclear polarisation, and magnetic-field effects. The theory combines spin-selective kinetics with singlet-triplet interconversion.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 13-14: Isotopes; implemented by `sys.isotopes={'1H','19F'}`.
- Lines 16-17: Chemical shifts; implemented by `inter.zeeman.scalar={0.0, 0.0}`.
- Lines 19-20: Chemical shift anisotropies (DFT); implemented by `inter.zeeman.eigs{1}=[0 0 0]`.
- Lines 25-26: Coordinates (DFT); implemented by `inter.coordinates={[0 0 0],[0 2.60 0]}`.
- Lines 28-29: J-coupling (expt); implemented by `inter.coupling.scalar=cell(2,2)`.
- Lines 32-33: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 39-40: Formalism and basis; implemented by `bas.formalism='sphten-liouv'`.
- Lines 43-44: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Get relaxation matrix; implemented by `R=relaxation(spin_system)`.
- Lines 50-51: Build relevant states; implemented by `U=unit_state(spin_system)`.
- Lines 58-59: Add pumping terms; implemented by `R=magpump(spin_system,R,Hz,1.3)`.
- Lines 62-63: Build active space projector; implemented by `P=[U Hz Fz -HzFz]`.
- Lines 65-66: Project the relaxation matrix into the active space; implemented by `disp('Kuprov''s matrix in Equation 2:'); disp(P'*R*P)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','19F'}`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0, 0.0}`.
- Lines 20: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[0 0 0]`.
- Lines 21: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[0 0 0]`.
- Lines 22: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[-47 -16 63]`.
- Lines 23: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[0 0 0]`.
- Lines 26: computes `inter.coordinates` using `inter.coordinates={[0 0 0],[0 2.60 0]}`.
- Lines 29: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2,2)`.
- Lines 30: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=50`.
- Lines 33: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 34: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 35: computes `inter.tau_c` using `inter.tau_c={110e-12}`.
- Lines 36: computes `inter.equilibrium` using `inter.equilibrium='IME'`.
- Lines 37: computes `inter.temperature` using `inter.temperature=298`.
- Lines 40: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 41: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 44: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- A simulation of the matrix in Equation 2 of IK's paper on chemically
- amplified NOEs (https://doi.org/10.1016/j.jmr.2004.01.011).
- Calculation time: seconds.
- Magnet field
- Isotopes
- Chemical shifts
- Chemical shift anisotropies (DFT)
- Coordinates (DFT)
- J-coupling (expt)
- Relaxation theory
- Formalism and basis
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `unit_state()`, `state()`, `magpump()`.
