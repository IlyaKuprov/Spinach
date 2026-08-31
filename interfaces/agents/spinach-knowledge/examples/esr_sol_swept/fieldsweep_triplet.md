# examples/esr_sol_swept/fieldsweep_triplet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_swept/fieldsweep_triplet.m`
- Signature: `fieldsweep_triplet()`
- Total lines: 64

## Purpose

Powder averaged X-band field-swept ESR spectrum of photo- generated pentacene triplet state. Calculation time: seconds.

## Physical / mathematical content

- Field-swept ESR examples. These files emphasise resonance-field finding, powder averaging, anisotropic g and hyperfine tensors, and intensity accumulation over orientation manifolds.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnet field (must be 1); implemented by `sys.magnet=1`.
- Lines 14-15: Triplet electron; implemented by `sys.isotopes={'E3'}`.
- Lines 17-18: Zeeman tensor, assumed isotropic; implemented by `inter.zeeman.matrix={diag([2.0 2.0 2.0])}`.
- Lines 20-21: ZFS, photo-excited pentacene triplet; implemented by `D=1360.1*1e6; E=-47.2*1e6`.
- Lines 24-25: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Experiment parameters; implemented by `parameters.spins={'E3'}`.
- Lines 43-44: Zeeman tensor into Hz/Tesla; implemented by `Z=-spin('E')*inter.zeeman.matrix{1,1}/(2*pi*2.00231930436256)`.
- Lines 46-52: Orientation-and field-dependent initial condition; implemented by `parameters.rho0=@(B,alp,bet,gam)zftrip(spin_system,euler2dcm(alp,bet,gam)* inter.coupling.matrix{1,1}* euler2dcm(alp,bet,gam)', [0.56 0.31 0.13], euler2dcm(alp,bet,gam)*…`.
- Lines 54-55: Run the simulation; implemented by `[spec,parameters]=fieldsweep(spin_system,parameters)`.
- Lines 57-58: Plotting; implemented by `kfigure(); plot(parameters.b_axis,spec)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=1`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'E3'}`.
- Lines 18: computes `inter.zeeman.matrix` using `inter.zeeman.matrix={diag([2.0 2.0 2.0])}`.
- Lines 21: computes `D` using `D=1360.1*1e6; E=-47.2*1e6`.
- Lines 22: computes `inter.coupling.matrix` using `inter.coupling.matrix={zfs2mat(D,E,0,0,0)}`.
- Lines 25: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `parameters.spins` using `parameters.spins={'E3'}`.
- Lines 34: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 35: computes `parameters.mw_freq` using `parameters.mw_freq=9e9`.
- Lines 36: computes `parameters.fwhm` using `parameters.fwhm=5e-4`.
- Lines 37: computes `parameters.int_tol` using `parameters.int_tol=0.1`.
- Lines 38: computes `parameters.tm_tol` using `parameters.tm_tol=0.1`.
- Lines 39: computes `parameters.window` using `parameters.window=[0.25 0.40]`.
- Lines 40: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 41: computes `parameters.rspt_order` using `parameters.rspt_order=Inf`.
- Lines 44: computes `Z` using `Z=-spin('E')*inter.zeeman.matrix{1,1}/(2*pi*2.00231930436256)`.

## Implementation structure

- Powder averaged X-band field-swept ESR spectrum of photo-
- generated pentacene triplet state.
- Calculation time: seconds.
- Magnet field (must be 1)
- Triplet electron
- Zeeman tensor, assumed isotropic
- ZFS, photo-excited pentacene triplet
- Basis set
- Spinach housekeeping
- Experiment parameters
- Zeeman tensor into Hz/Tesla
- Orientation-and field-dependent initial condition

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `zfs2mat()`, `create()`, `basis()`, `spin()`, `zftrip()`, `euler2dcm()`, `fieldsweep()`, `kfigure()`, `kxlabel()`, `kylabel()`.
