# examples/dnp_liq/odnp_liquid_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/dnp_liq/odnp_liquid_2.m`
- Signature: `odnp_liquid_2()`
- Total lines: 80

## Purpose

Overhauser type DNP in liquid phase at room temperature, after a perfect inversion pulse on the electron ESR signal. The simulation uses Redfield theory to account for the dipolar cross-relaxation. Calculation time: seconds

## Physical / mathematical content

- Liquid-state DNP examples. The main ingredients are electron-nuclear cross-relaxation, scalar or dipolar contact mechanisms, motional spectral densities, and field/frequency dependence of polarisation transfer.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Spin system; implemented by `sys.magnet=3.4`.
- Lines 15-18: Zeeman interactions; implemented by `inter.zeeman.matrix={[5 0 0; 0 5 0; 0 0 5] [5 0 0; 0 5 0; 0 0 5] [2.0023 0 0; 0 2.0025 0; 0 0 2.0027]}`.
- Lines 20-21: Coordinates (Angstrom); implemented by `inter.coordinates={[0.0 0.0 0.0]`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 29-30: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Isotropic thermal equilibrium; implemented by `rho_eq=equilibrium(spin_system)`.
- Lines 43-44: Electron control operator; implemented by `Lx=operator(spin_system,'Lx','E')`.
- Lines 46-47: Inversion pulse on the electron; implemented by `rho0=step(spin_system,Lx,rho_eq,pi)`.
- Lines 49-50: Experiment paramaters; implemented by `parameters.spins={'E'}`.
- Lines 62-63: Simulation; implemented by `answer=liquid(spin_system,@dnp_time_dep,parameters,'esr')`.
- Lines 65-66: Plotting; implemented by `kfigure(); x_axis=linspace(0,1000,1001)`.

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
- Lines 41: computes `rho_eq` using `rho_eq=equilibrium(spin_system)`.
- Lines 44: computes `Lx` using `Lx=operator(spin_system,'Lx','E')`.
- Lines 47: computes `rho0` using `rho0=step(spin_system,Lx,rho_eq,pi)`.
- Lines 50: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 51: computes `parameters.rho0` using `parameters.rho0=rho0`.
- Lines 52-54: computes `parameters.coil` using `parameters.coil=[state(spin_system,{'Lz'},{1}) state(spin_system,{'Lz'},{2}) state(spin_system,{'Lz'},{3})]`.

## Implementation structure

- Overhauser type DNP in liquid phase at room temperature, after a perfect
- inversion pulse on the electron ESR signal. The simulation uses Redfield
- theory to account for the dipolar cross-relaxation.
- Calculation time: seconds
- Spin system
- Zeeman interactions
- Coordinates (Angstrom)
- Basis set
- Relaxation theory
- Spinach housekeeping
- Isotropic thermal equilibrium
- Electron control operator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `equilibrium()`, `operator()`, `step()`, `state()`, `liquid()`, `kfigure()`, `subplot()`, `answer()`, `kxlabel()`, `kylabel()`, `klegend()`.
