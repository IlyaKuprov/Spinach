# examples/imaging/press_3d_example.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/press_3d_example.m`
- Signature: `press_3d_example()`
- Total lines: 79

## Purpose

PRESS excitation profile in three dimensions with tilted gradient system. Change the frequency under Pulse Parame- ters to move the hot spot through the sample. Simulation time: hours, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnetic induction; implemented by `sys.magnet=3.0`.
- Lines 15-16: Spin systems; implemented by `sys.isotopes={'1H','1H'}`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Disable path tracing; implemented by `sys.disable={'pt'}`.
- Lines 28-29: This is here essential; implemented by `sys.enable={'polyadic'}`.
- Lines 31-32: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 35-36: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 48-49: Pulse parameters; implemented by `parameters.rf_phi={pi/2 pi/2 pi/2}`.
- Lines 55-56: Sample geometry; implemented by `parameters.dims=[0.30 0.25 0.27]`.
- Lines 61-62: Relaxation phantom; implemented by `parameters.rlx_ph={zeros(parameters.npts)}`.
- Lines 65-66: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts)}`.
- Lines 71-72: Run voxel selection diagnostics; implemented by `phan=imaging(spin_system,@press_voxel_3d,parameters)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=3.0`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.5,3.7}`.
- Lines 18: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2)`.
- Lines 19: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=10`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 29: computes `sys.enable` using `sys.enable={'polyadic'}`.
- Lines 32: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 36: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 37: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 38: computes `parameters.offset` using `parameters.offset=0`.
- Lines 39: computes `parameters.sweep` using `parameters.sweep=1000`.
- Lines 40: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 41: computes `parameters.npoints` using `parameters.npoints=128`.
- Lines 42: computes `parameters.invert_axis` using `parameters.invert_axis=1`.
- Lines 43: computes `parameters.ss_grad_amp` using `parameters.ss_grad_amp=[25e-3 25e-3 25e-3]`.

## Implementation structure

- PRESS excitation profile in three dimensions with tilted
- gradient system. Change the frequency under Pulse Parame-
- ters to move the hot spot through the sample.
- Simulation time: hours, faster with a Tesla V100 GPU.
- Magnetic induction
- Spin systems
- Basis set
- Disable path tracing
- This is here essential
- Spinach housekeeping
- Sequence parameters
- Pulse parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `state()`, `imaging()`, `volplot()`, `kxlabel()`, `kylabel()`, `kzlabel()`.
