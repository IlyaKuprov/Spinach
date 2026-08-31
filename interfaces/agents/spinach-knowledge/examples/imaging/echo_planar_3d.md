# examples/imaging/echo_planar_3d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/echo_planar_3d.m`
- Signature: `echo_planar_3d()`
- Total lines: 125

## Purpose

Slice selection in 3D followed by three-dimensional echo planar imaging sequence. Simulation time: hours, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Isotopes; implemented by `sys.isotopes={'1H'}`.
- Lines 14-15: Magnetic induction; implemented by `sys.magnet=5.9`.
- Lines 17-18: Chemical shifts; implemented by `inter.zeeman.scalar={0.0}`.
- Lines 20-21: Relaxation model; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 27-28: Disable path tracing; implemented by `sys.disable={'pt'}`.
- Lines 30-31: This needs a GPU; implemented by `sys.enable={'greedy'}`.
- Lines 33-34: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 41-42: Pulse phase; implemented by `parameters.rf_phi=pi/2`.
- Lines 44-45: Number of steps in the pulse; implemented by `pulse_nsteps=50`.
- Lines 47-48: Overall pulse duration; implemented by `pulse_time=2.0e-4`.
- Lines 50-51: Slice selection pulse frequency table; implemented by `parameters.rf_frq_list=-5e3*ones(1,pulse_nsteps)`.
- Lines 53-54: Slice selection pulse amplitude table; implemented by `parameters.rf_amp_list=2*pi*7500*pulse_shape('gaussian',pulse_nsteps)`.
- Lines 56-57: Slice selection pulse duration table; implemented by `parameters.rf_dur_list=(pulse_time/pulse_nsteps)*ones(1,pulse_nsteps)`.
- Lines 59-60: Sequence parameters; implemented by `parameters.image_size=[201 201]`.
- Lines 69-70: Phantom library call; implemented by `[R1_Ph,R2_Ph,PD_Ph,dims,npts]=phantoms('brain-medres')`.
- Lines 72-73: Sample settings; implemented by `parameters.dims=dims`.
- Lines 79-80: Relaxation phantom; implemented by `[R1,R2]=rlx_t1_t2(spin_system)`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 15: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 21: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 22: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 23: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 24: computes `inter.r1_rates` using `inter.r1_rates={1}`.
- Lines 25: computes `inter.r2_rates` using `inter.r2_rates={1}`.
- Lines 28: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 31: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 34: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 35: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 38: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 42: computes `parameters.rf_phi` using `parameters.rf_phi=pi/2`.
- Lines 45: computes `pulse_nsteps` using `pulse_nsteps=50`.
- Lines 48: computes `pulse_time` using `pulse_time=2.0e-4`.
- Lines 51: computes `parameters.rf_frq_list` using `parameters.rf_frq_list=-5e3*ones(1,pulse_nsteps)`.
- Lines 54: computes `parameters.rf_amp_list` using `parameters.rf_amp_list=2*pi*7500*pulse_shape('gaussian',pulse_nsteps)`.

## Implementation structure

- Slice selection in 3D followed by three-dimensional echo
- planar imaging sequence.
- Simulation time: hours, faster with a Tesla V100 GPU.
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation model
- Disable path tracing
- This needs a GPU
- Basis set
- Spinach housekeeping
- Pulse phase

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `pulse_shape()`, `phantoms()`, `rlx_t1_t2()`, `state()`, `kfigure()`, `dims()`, `volplot()`, `ktitle()`, `imaging()`, `scale_figure()`, `subplot()`, `mri_2d_plot()`, `apodisation()`, `fftshift()`.
