# examples/esr_sol_swept/temperature_gd_dota.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_swept/temperature_gd_dota.m`
- Signature: `temperature_gd_dota()`
- Total lines: 67

## Purpose

Powder averaged W-band field-swept ESR spectrum of Gd(III) DOTA complex. Exact diagonalisation is used and a tempera- ture dependence plot is produced. Calculation time: seconds.

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
- Lines 26-27: Get figure going; implemented by `kfigure(); scale_figure([2.0 1.5])`.
- Lines 29-30: Temperatures; implemented by `T=[100 10 1 0.1]`.
- Lines 32-33: Loop over temperatures; implemented by `for n=1:4`.
- Lines 35-36: Set the temperature; implemented by `inter.temperature=T(n)`.
- Lines 38-39: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Experiment parameters; implemented by `parameters.spins={'E8'}`.
- Lines 53-54: Run the simulation; implemented by `[spec,parameters]=fieldsweep(spin_system,parameters)`.
- Lines 56-57: Plotting; implemented by `subplot(2,2,n); plot(parameters.b_axis,spec)`.

### Control flow inferred from the code

- Line 33: `for` loop over `n=1:4`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'E8'}`.
- Lines 15: computes `sys.magnet` using `sys.magnet=1`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.9918}`.
- Lines 19: computes `inter.coupling.eigs{1,1}` using `inter.coupling.eigs{1,1}=[0.57e9 0.57e9 -2*0.57e9]/3`.
- Lines 20: computes `inter.coupling.euler{1,1}` using `inter.coupling.euler{1,1}=[0 0 0]`.
- Lines 23: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `T` using `T=[100 10 1 0.1]`.
- Lines 36: computes `inter.temperature` using `inter.temperature=T(n)`.
- Lines 39: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 43: computes `parameters.spins` using `parameters.spins={'E8'}`.
- Lines 44: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 45: computes `parameters.mw_freq` using `parameters.mw_freq=90e9`.
- Lines 46: computes `parameters.fwhm` using `parameters.fwhm=2e-4`.
- Lines 47: computes `parameters.int_tol` using `parameters.int_tol=1.0`.
- Lines 48: computes `parameters.tm_tol` using `parameters.tm_tol=0.1`.
- Lines 49: computes `parameters.window` using `parameters.window=[3.05 3.4]`.
- Lines 50: computes `parameters.npoints` using `parameters.npoints=4096`.

## Implementation structure

- Powder averaged W-band field-swept ESR spectrum of Gd(III)
- DOTA complex. Exact diagonalisation is used and a tempera-
- ture dependence plot is produced.
- Calculation time: seconds.
- Isotopes
- Magnet field (must be 1)
- Properties
- Basis set
- Get figure going
- Temperatures
- Loop over temperatures
- Set the temperature

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `kfigure()`, `scale_figure()`, `create()`, `basis()`, `fieldsweep()`, `subplot()`, `kxlabel()`, `kylabel()`, `ktitle()`, `num2str()`.
