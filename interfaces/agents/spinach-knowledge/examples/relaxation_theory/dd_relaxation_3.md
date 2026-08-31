# examples/relaxation_theory/dd_relaxation_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/relaxation_theory/dd_relaxation_3.m`
- Signature: `dd_relaxation_3()`
- Total lines: 76

## Purpose

Extreme narrowing limit case comparison between the dipolar relaxation rates in proton-proton and proton-deuterium system. The rate must sca- le with the square of the magnetogyric ratio and with S(S+1), where S is the quantum number of the partner spin. Calculation time: seconds

## Physical / mathematical content

- Relaxation-theory examples. The mathematical backbone is Bloch-Redfield-Wangsness or stochastic Liouville theory, spectral densities, cross-correlation terms, motional models, and extraction of longitudinal/transverse decay behaviour from superoperators.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-15: % System with two protons; implemented by `sys.magnet=14.1`.
- Lines 14-15: System specification; implemented by `sys.magnet=14.1`.
- Lines 20-21: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Relaxation superoperator; implemented by `R=relaxation(spin_system)`.
- Lines 37-38: Longitudinal relaxation rate; implemented by `rho=state(spin_system,{'Lz'},{1})`.
- Lines 41-44: % System with a proton and a deuteron; implemented by `sys.magnet=14.1`.
- Lines 43-44: System specification -two protons; implemented by `sys.magnet=14.1`.
- Lines 70-71: % Comparison; implemented by `disp(['R1 rate ratio, Spinach: ' num2str(R1pp/R1pd)])`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 17: computes `inter.coordinates` using `inter.coordinates={[0.0 0.0 0.0]`.
- Lines 21: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 22: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 23: computes `inter.rlx_keep` using `inter.rlx_keep='labframe'`.
- Lines 24: computes `inter.tau_c` using `inter.tau_c={1e-12}`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `R` using `R=relaxation(spin_system)`.
- Lines 38: computes `rho` using `rho=state(spin_system,{'Lz'},{1})`.
- Lines 72: computes `ratio` using `ratio=(((1/2)*(1/2+1))/(1*(1+1)))*(spin('1H')/spin('2H'))^2`.

## Implementation structure

- Extreme narrowing limit case comparison between the dipolar relaxation
- rates in proton-proton and proton-deuterium system. The rate must sca-
- le with the square of the magnetogyric ratio and with S(S+1), where S
- is the quantum number of the partner spin.
- Calculation time: seconds
- % System with two protons
- System specification
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Relaxation superoperator
- Longitudinal relaxation rate

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `state()`, `num2str()`, `spin()`.
