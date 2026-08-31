# examples/relaxation_theory/dd_relaxation_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/dd_relaxation_2.m`
- Signature: `dd_relaxation_2()`
- Total lines: 42

## Purpose

Complete Bloch-Redfield-Wangsness relaxation superoperator in a system with dipolar coupling between spins. The dipolar couplings are computed from Cartesian coordinates of the spins. The result should not depend on the choice of the rotation angles below. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: System specification; implemented by `sys.magnet=14.1`.
- Lines 17-18: Randomly rotated set of coordinates; implemented by `R=euler2dcm([pi/3 pi/4 pi/5])`.
- Lines 23-24: Relaxation theory parameters; implemented by `sys.tols.rlx_integration=1e-5`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Relaxation superoperator; implemented by `disp(full(relaxation(spin_system)))`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.0,2.0,3.0}`.
- Lines 18: computes `R` using `R=euler2dcm([pi/3 pi/4 pi/5])`.
- Lines 19: computes `inter.coordinates` using `inter.coordinates={[-0.2230 0.7893 -0.5721]*R`.
- Lines 24: computes `sys.tols.rlx_integration` using `sys.tols.rlx_integration=1e-5`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 26: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 27: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 28: computes `inter.tau_c` using `inter.tau_c={1e-9}`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- Complete Bloch-Redfield-Wangsness relaxation superoperator in a system
- with dipolar coupling between spins. The dipolar couplings are computed
- from Cartesian coordinates of the spins. The result should not depend
- on the choice of the rotation angles below.
- Calculation time: seconds
- System specification
- Randomly rotated set of coordinates
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Relaxation superoperator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `euler2dcm()`, `create()`, `basis()`, `relaxation()`.
