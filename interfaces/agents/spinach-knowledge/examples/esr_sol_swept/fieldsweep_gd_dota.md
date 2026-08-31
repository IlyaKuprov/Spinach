# examples/esr_sol_swept/fieldsweep_gd_dota.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_swept/fieldsweep_gd_dota.m`
- Signature: `fieldsweep_gd_dota()`
- Total lines: 52

## Purpose

Powder averaged W-band field-swept ESR spectrum of Gd(III) DOTA complex. Exact diagonalisation is used. Calculation time: seconds.

## Physical / mathematical content

- Field-swept ESR examples. These files emphasise resonance-field finding, powder averaging, anisotropic g and hyperfine tensors, and intensity accumulation over orientation manifolds.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Isotopes; implemented by `sys.isotopes={'E8'}`.
- Lines 14-15: Magnet field (must be 1); implemented by `sys.magnet=1`.
- Lines 17-18: Properties; implemented by `inter.zeeman.scalar={1.9918}`.
- Lines 22-23: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 26-27: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Experiment parameters; implemented by `parameters.spins={'E8'}`.
- Lines 41-42: Run the simulation in the high-T approximation; implemented by `parameters.rho0=-state(spin_system,'Lz','E8')`.
- Lines 45-46: Plotting; implemented by `kfigure(); plot(parameters.b_axis,spec)`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'E8'}`.
- Lines 15: computes `sys.magnet` using `sys.magnet=1`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.9918}`.
- Lines 19: computes `inter.coupling.eigs{1,1}` using `inter.coupling.eigs{1,1}=[0.57e9 0.57e9 -2*0.57e9]/3`.
- Lines 20: computes `inter.coupling.euler{1,1}` using `inter.coupling.euler{1,1}=[0 0 0]`.
- Lines 23: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `parameters.spins` using `parameters.spins={'E8'}`.
- Lines 32: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 33: computes `parameters.mw_freq` using `parameters.mw_freq=90e9`.
- Lines 34: computes `parameters.fwhm` using `parameters.fwhm=2e-4`.
- Lines 35: computes `parameters.int_tol` using `parameters.int_tol=10.0`.
- Lines 36: computes `parameters.tm_tol` using `parameters.tm_tol=0.1`.
- Lines 37: computes `parameters.window` using `parameters.window=[3.05 3.4]`.
- Lines 38: computes `parameters.npoints` using `parameters.npoints=4096`.
- Lines 39: computes `parameters.rspt_order` using `parameters.rspt_order=Inf`.
- Lines 42: computes `parameters.rho0` using `parameters.rho0=-state(spin_system,'Lz','E8')`.

## Implementation structure

- Powder averaged W-band field-swept ESR spectrum of Gd(III)
- DOTA complex. Exact diagonalisation is used.
- Calculation time: seconds.
- Isotopes
- Magnet field (must be 1)
- Properties
- Basis set
- Spinach housekeeping
- Experiment parameters
- Run the simulation in the high-T approximation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `fieldsweep()`, `kfigure()`, `kxlabel()`, `kylabel()`.
