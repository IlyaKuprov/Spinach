# examples/dnp_liq/ccdnp/rates_si_sys_b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/ccdnp/rates_si_sys_b.m`
- Signature: `rates_si_sys_b()`
- Total lines: 100

## Purpose

Self-and cross-relaxation rates in cross-correlated DNP, considering a system with two electrons connected by exchange coupling, both cou- pled to a nucleus by dipolar couplings. Further particulars in: Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 17-18: Spin system; implemented by `sys.isotopes={'1H','E','E'}`.
- Lines 20-21: Zeeman interactions; implemented by `inter.zeeman.eigs{1,1}=[0 10 20]`.
- Lines 28-29: Exchange coupling; implemented by `inter.coupling.scalar=cell(3,3)`.
- Lines 32-33: coordinates; implemented by `inter.coordinates={[ 0.000 0.000 0.000]`.
- Lines 37-38: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 41-42: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 49-50: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 53-54: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 56-57: Rate printouts; implemented by `disp(' ')`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'1H','E','E'}`.
- Lines 21: computes `inter.zeeman.eigs{1,1}` using `inter.zeeman.eigs{1,1}=[0 10 20]`.
- Lines 22: computes `inter.zeeman.euler{1,1}` using `inter.zeeman.euler{1,1}=[0 0 0]`.
- Lines 23: computes `inter.zeeman.eigs{1,2}` using `inter.zeeman.eigs{1,2}=[1.977800 1.977600 1.977600]`.
- Lines 24: computes `inter.zeeman.euler{1,2}` using `inter.zeeman.euler{1,2}=[-0.180 0.017 0.194]`.
- Lines 25: computes `inter.zeeman.eigs{1,3}` using `inter.zeeman.eigs{1,3}=[2.006800 2.003800 2.003800]`.
- Lines 26: computes `inter.zeeman.euler{1,3}` using `inter.zeeman.euler{1,3}=[0.632 0.783 1.086]`.
- Lines 29: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(3,3)`.
- Lines 30: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=5e6`.
- Lines 33: computes `inter.coordinates` using `inter.coordinates={[ 0.000 0.000 0.000]`.
- Lines 38: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 39: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 42: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 43: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 44: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 45: computes `inter.temperature` using `inter.temperature=298`.
- Lines 46: computes `inter.tau_c` using `inter.tau_c={100e-12}`.

## Implementation structure

- Self-and cross-relaxation rates in cross-correlated DNP, considering
- a system with two electrons connected by exchange coupling, both cou-
- pled to a nucleus by dipolar couplings. Further particulars in:
- Calculation time: seconds
- Magnet field
- Spin system
- Zeeman interactions
- Exchange coupling
- coordinates
- Basis set
- Relaxation theory
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `state()`, `num2str()`.
