# examples/extremes/ph_enc_3d_highres.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/extremes/ph_enc_3d_highres.m`
- Signature: `ph_enc_3d_highres()`
- Total lines: 124

## Purpose

Slice selection in 3D followed by phase-encoded imaging of the resulting slice. This simulation fills up a sys- tem with eight H200 GPUs and 4 TB of RAM. Simulation time: you hope and pray this even starts; if it does, then hours.

## Physical / mathematical content

- Extreme-regime examples. These scripts exercise Spinach in unusually large, stiff, high-field, low-field, or otherwise numerically demanding regimes where approximations, conditioning, and basis-size control are central.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Isotopes; implemented by `sys.isotopes={'1H'}`.
- Lines 15-16: Magnetic induction; implemented by `sys.magnet=5.9`.
- Lines 18-19: Chemical shifts; implemented by `inter.zeeman.scalar={0.0}`.
- Lines 21-22: Relaxation theory; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 28-29: Disable path tracing; implemented by `sys.disable={'pt'}`.
- Lines 31-32: This needs a GPU; implemented by `sys.enable={'greedy','polyadic'}`.
- Lines 34-35: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 38-39: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Gat phantom from library; implemented by `[R1_Ph,R2_Ph,PD_Ph,dims,npts]=phantoms('brain-highres')`.
- Lines 45-46: Sample settings; implemented by `parameters.dims=dims`.
- Lines 52-53: Relaxation operators and phantoms; implemented by `[R1,R2]=rlx_t1_t2(spin_system)`.
- Lines 57-58: Initial and detection state phantoms; implemented by `parameters.rho0_ph={PD_Ph}`.
- Lines 63-64: Diffusion and flow; implemented by `parameters.u=zeros(parameters.npts)`.
- Lines 69-70: Pulse phase; implemented by `parameters.rf_phi=pi/2`.
- Lines 72-73: Number of steps in the pulse; implemented by `pulse_nsteps=50`.
- Lines 75-76: Overall pulse duration; implemented by `pulse_time=2.0e-4`.
- Lines 78-79: Slice selection pulse frequency table; implemented by `parameters.rf_frq_list=-5e3*ones(1,pulse_nsteps)`.
- Lines 81-82: Slice selection pulse amplitude table; implemented by `parameters.rf_amp_list=2*pi*7500*pulse_shape('gaussian',pulse_nsteps)`.

### Key state/data transformations

- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 16: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 19: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 22: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 23: computes `inter.r1_rates` using `inter.r1_rates={1.0}`.
- Lines 24: computes `inter.r2_rates` using `inter.r2_rates={1.0}`.
- Lines 25: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 29: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 32: computes `sys.enable` using `sys.enable={'greedy','polyadic'}`.
- Lines 35: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 36: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 39: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 43: computes `[R1_Ph,R2_Ph,PD_Ph,dims,npts]` using `[R1_Ph,R2_Ph,PD_Ph,dims,npts]=phantoms('brain-highres')`.
- Lines 46: computes `parameters.dims` using `parameters.dims=dims`.
- Lines 47: computes `parameters.npts` using `parameters.npts=npts`.
- Lines 48: computes `parameters.deriv` using `parameters.deriv={'period',3}`.
- Lines 49: computes `parameters.spins` using `parameters.spins={'1H'}`.

## Implementation structure

- Slice selection in 3D followed by phase-encoded imaging
- of the resulting slice. This simulation fills up a sys-
- tem with eight H200 GPUs and 4 TB of RAM.
- Simulation time: you hope and pray this even starts;
- if it does, then hours.
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation theory
- Disable path tracing
- This needs a GPU
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `phantoms()`, `rlx_t1_t2()`, `state()`, `pulse_shape()`, `kfigure()`, `dims()`, `volplot()`, `ktitle()`, `imaging()`, `scale_figure()`, `subplot()`, `mri_2d_plot()`, `apodisation()`, `fftshift()`.
