# examples/esr_sol_swept/fieldsweep_nitroxide.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_swept/fieldsweep_nitroxide.m`
- Signature: `fieldsweep_nitroxide()`
- Total lines: 56

## Purpose

Field swept EPR spectrum of nitroxide, computed by finding resonance fields and transition moments. Calculation time: seconds.

## Physical / mathematical content

- Field-swept ESR examples. These files emphasise resonance-field finding, powder averaging, anisotropic g and hyperfine tensors, and intensity accumulation over orientation manifolds.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Isotopes; implemented by `sys.isotopes={'E','14N'}`.
- Lines 13-14: Magnet field (must be 1); implemented by `sys.magnet=1`.
- Lines 16-17: Interactions; implemented by `inter.zeeman.matrix=cell(1,2)`.
- Lines 26-27: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Experiment parameters; implemented by `parameters.spins={'E'}`.
- Lines 45-46: Run the simulation in the high-T approximation; implemented by `parameters.rho0=-state(spin_system,'Lz','E')`.
- Lines 49-50: Plotting; implemented by `kfigure(); plot(parameters.b_axis,spec)`.

### Key state/data transformations

- Lines 11: computes `sys.isotopes` using `sys.isotopes={'E','14N'}`.
- Lines 14: computes `sys.magnet` using `sys.magnet=1`.
- Lines 17: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,2)`.
- Lines 18: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=[2.01045 0.00000 0.00000`.
- Lines 21: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(2,2)`.
- Lines 22: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=[1.2356 0.0000 0.6322`.
- Lines 27: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 36: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 37: computes `parameters.mw_freq` using `parameters.mw_freq=9e9`.
- Lines 38: computes `parameters.fwhm` using `parameters.fwhm=1e-5`.
- Lines 39: computes `parameters.int_tol` using `parameters.int_tol=10.0`.
- Lines 40: computes `parameters.tm_tol` using `parameters.tm_tol=0.1`.
- Lines 41: computes `parameters.window` using `parameters.window=[0.316 0.326]`.
- Lines 42: computes `parameters.npoints` using `parameters.npoints=1024`.
- Lines 43: computes `parameters.rspt_order` using `parameters.rspt_order=Inf`.

## Implementation structure

- Field swept EPR spectrum of nitroxide, computed by finding
- resonance fields and transition moments.
- Calculation time: seconds.
- Isotopes
- Magnet field (must be 1)
- Interactions
- Basis set
- Spinach housekeeping
- Experiment parameters
- Run the simulation in the high-T approximation
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `fieldsweep()`, `kfigure()`, `kxlabel()`, `kylabel()`.
