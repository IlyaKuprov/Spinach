# examples/dnp_liq/odnp_liquid_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/odnp_liquid_1.m`
- Signature: `odnp_liquid_1()`
- Total lines: 69

## Purpose

Overhauser type DNP in liquid phase at room temperature, using a continu- ous on-resonance CW irradiation of the electron ESR signal. The simulati- on uses Redfield theory to account for the dipolar cross-relaxation. Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Spin system; implemented by `sys.magnet=3.4`.
- Lines 15-18: Zeeman interactions; implemented by `inter.zeeman.matrix={[5 0 0; 0 5 0; 0 0 5] [5 0 0; 0 5 0; 0 0 5] [2.0023 0 0; 0 2.0025 0; 0 0 2.0027]}`.
- Lines 20-21: Coordinates (Angstrom); implemented by `inter.coordinates={[0.0 0.0 0.0]`.
- Lines 25-26: Complete basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 29-30: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Experiment paramaters; implemented by `parameters.spins={'E'}`.
- Lines 53-54: Simulation; implemented by `answer=liquid(spin_system,@dnp_time_dep,parameters,'esr')`.
- Lines 56-57: Plotting; implemented by `kfigure(); x_axis=linspace(0,1000,1001)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=3.4`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H','1H','E'}`.
- Lines 16-18: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={[5 0 0; 0 5 0; 0 0 5] [5 0 0; 0 5 0; 0 0 5] [2.0023 0 0; 0 2.0025 0; 0 0 2.0027]}`.
- Lines 21: computes `inter.coordinates` using `inter.coordinates={[0.0 0.0 0.0]`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 31: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 32: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 33: computes `inter.temperature` using `inter.temperature=298`.
- Lines 34: computes `inter.tau_c` using `inter.tau_c={10e-12}`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 42: computes `parameters.needs` using `parameters.needs={'rho_eq'}`.
- Lines 43-45: computes `parameters.coil` using `parameters.coil=[state(spin_system,{'Lz'},{1}) state(spin_system,{'Lz'},{2}) state(spin_system,{'Lz'},{3})]`.
- Lines 46: computes `parameters.mw_pwr` using `parameters.mw_pwr=2*pi*1e6`.
- Lines 47: computes `parameters.mw_off` using `parameters.mw_off=0`.
- Lines 48: computes `parameters.mw_oper` using `parameters.mw_oper=operator(spin_system,'Lx','E')`.

## Implementation structure

- Overhauser type DNP in liquid phase at room temperature, using a continu-
- ous on-resonance CW irradiation of the electron ESR signal. The simulati-
- on uses Redfield theory to account for the dipolar cross-relaxation.
- Calculation time: seconds
- Spin system
- Zeeman interactions
- Coordinates (Angstrom)
- Complete basis set
- Relaxation theory
- Spinach housekeeping
- Experiment paramaters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `liquid()`, `kfigure()`, `subplot()`, `answer()`, `kxlabel()`, `kylabel()`, `klegend()`.
