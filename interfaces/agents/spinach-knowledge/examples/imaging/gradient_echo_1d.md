# examples/imaging/gradient_echo_1d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/gradient_echo_1d.m`
- Signature: `gradient_echo_1d()`
- Total lines: 64

## Purpose

A gradient echo experiment in the presence of diffusion and flow. Calculation time: seconds. Ahmed Allami Ilya Kuprov

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Isotopes; implemented by `sys.isotopes={'1H'}`.
- Lines 14-15: Magnetic induction; implemented by `sys.magnet=5.9`.
- Lines 17-18: Chemical shifts; implemented by `inter.zeeman.scalar={1.0}`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 35-36: Sample geometry; implemented by `parameters.dims=0.30`.
- Lines 40-41: Relaxation phantom; implemented by `parameters.rlx_ph={}`.
- Lines 44-45: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 50-51: Diffusion and flow; implemented by `parameters.u=ones(parameters.npts,1)`.
- Lines 54-55: Run the simulation; implemented by `echo=imaging(spin_system,@grad_echo,parameters)`.
- Lines 57-58: Plotting; implemented by `grad_duration=parameters.g_step_dur*parameters.g_n_steps`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 15: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.0}`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 29: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 30: computes `parameters.offset` using `parameters.offset=0.0`.
- Lines 31: computes `parameters.g_amp` using `parameters.g_amp=5e-6`.
- Lines 32: computes `parameters.g_step_dur` using `parameters.g_step_dur=2e-4`.
- Lines 33: computes `parameters.g_n_steps` using `parameters.g_n_steps=200`.
- Lines 36: computes `parameters.dims` using `parameters.dims=0.30`.
- Lines 37: computes `parameters.npts` using `parameters.npts=100`.
- Lines 38: computes `parameters.deriv` using `parameters.deriv={'period',3}`.
- Lines 41: computes `parameters.rlx_ph` using `parameters.rlx_ph={}`.
- Lines 42: computes `parameters.rlx_op` using `parameters.rlx_op={}`.
- Lines 45: computes `parameters.rho0_ph` using `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 46: computes `parameters.rho0_st` using `parameters.rho0_st={state(spin_system,'Lz','1H')}`.

## Implementation structure

- A gradient echo experiment in the presence
- of diffusion and flow.
- Calculation time: seconds.
- Ahmed Allami
- Ilya Kuprov
- Isotopes
- Magnetic induction
- Chemical shifts
- Basis set
- Spinach housekeeping
- Sequence parameters
- Sample geometry

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `imaging()`, `kfigure()`, `kxlabel()`, `kylabel()`.
