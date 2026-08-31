# examples/relaxation_theory/scalar_relaxation_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/scalar_relaxation_1.m`
- Signature: `scalar_relaxation_1()`
- Total lines: 43

## Purpose

Redfield superoperator for the scalar relaxation of the first kind in a two-proton system with a noisy J-coupling. This si- tuation occurs in aziridines, where the slow nitrogen inversi- on jitters scalar couplings on a millisecond time scale. Set to demonstrate the effect described in: Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: System specification; implemented by `sys.magnet=11.75`.
- Lines 22-23: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 26-27: Relaxation superoperator; implemented by `inter.relaxation={'SRFK'}`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Show a spy plot of R; implemented by `kfigure(); spy(relaxation(spin_system))`.

### Key state/data transformations

- Lines 18: computes `sys.magnet` using `sys.magnet=11.75`.
- Lines 19: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 20: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0 2.0}`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `inter.relaxation` using `inter.relaxation={'SRFK'}`.
- Lines 28: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 29: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 30: computes `inter.srfk_tau_c` using `inter.srfk_tau_c={[1.0 1e-3]}`.
- Lines 31: computes `inter.srfk_mdepth` using `inter.srfk_mdepth=cell(2)`.
- Lines 32: computes `inter.srfk_mdepth{1,2}` using `inter.srfk_mdepth{1,2}=15.0`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- Redfield superoperator for the scalar relaxation of the first
- kind in a two-proton system with a noisy J-coupling. This si-
- tuation occurs in aziridines, where the slow nitrogen inversi-
- on jitters scalar couplings on a millisecond time scale. Set
- to demonstrate the effect described in:
- Calculation time: seconds
- System specification
- Basis set
- Relaxation superoperator
- Spinach housekeeping
- Show a spy plot of R

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `kfigure()`, `spy()`, `relaxation()`, `ktitle()`.
