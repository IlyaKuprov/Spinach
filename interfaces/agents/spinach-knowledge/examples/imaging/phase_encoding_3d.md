# examples/imaging/phase_encoding_3d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/phase_encoding_3d.m`
- Signature: `phase_encoding_3d()`
- Total lines: 123

## Purpose

Slice selection in 3D followed by phase-encoded imaging of the resulting slice. Simulation time: minutes, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Isotopes; implemented by `sys.isotopes={'1H'}`.
- Lines 14-15: Magnetic induction; implemented by `sys.magnet=5.9`.
- Lines 17-18: Chemical shifts; implemented by `inter.zeeman.scalar={0.0}`.
- Lines 20-21: Relaxation theory; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 27-28: Disable path tracing; implemented by `sys.disable={'pt','krylov'}`.
- Lines 30-31: This needs a GPU; implemented by `sys.enable={'greedy'}`.
- Lines 33-34: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 41-42: Gat phantom from library; implemented by `[R1_Ph,R2_Ph,PD_Ph,dims,npts]=phantoms('brain-medres')`.
- Lines 44-45: Sample settings; implemented by `parameters.dims=dims`.
- Lines 51-52: Relaxation operators and phantoms; implemented by `[R1,R2]=rlx_t1_t2(spin_system)`.
- Lines 56-57: Initial and detection state phantoms; implemented by `parameters.rho0_ph={PD_Ph}`.
- Lines 62-63: Diffusion and flow; implemented by `parameters.u=zeros(parameters.npts)`.
- Lines 68-69: Pulse phase; implemented by `parameters.rf_phi=pi/2`.
- Lines 71-72: Number of steps in the pulse; implemented by `pulse_nsteps=50`.
- Lines 74-75: Overall pulse duration; implemented by `pulse_time=2.0e-4`.
- Lines 77-78: Slice selection pulse frequency table; implemented by `parameters.rf_frq_list=-5e3*ones(1,pulse_nsteps)`.
- Lines 80-81: Slice selection pulse amplitude table; implemented by `parameters.rf_amp_list=2*pi*7500*pulse_shape('gaussian',pulse_nsteps)`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 15: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 21: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 22: computes `inter.r1_rates` using `inter.r1_rates={1.0}`.
- Lines 23: computes `inter.r2_rates` using `inter.r2_rates={1.0}`.
- Lines 24: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 25: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 28: computes `sys.disable` using `sys.disable={'pt','krylov'}`.
- Lines 31: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 34: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 35: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 38: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 42: computes `[R1_Ph,R2_Ph,PD_Ph,dims,npts]` using `[R1_Ph,R2_Ph,PD_Ph,dims,npts]=phantoms('brain-medres')`.
- Lines 45: computes `parameters.dims` using `parameters.dims=dims`.
- Lines 46: computes `parameters.npts` using `parameters.npts=npts`.
- Lines 47: computes `parameters.deriv` using `parameters.deriv={'period',3}`.
- Lines 48: computes `parameters.spins` using `parameters.spins={'1H'}`.

## Implementation structure

- Slice selection in 3D followed by phase-encoded imaging
- of the resulting slice.
- Simulation time: minutes, faster with a Tesla V100 GPU.
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation theory
- Disable path tracing
- This needs a GPU
- Basis set
- Spinach housekeeping
- Gat phantom from library

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `phantoms()`, `rlx_t1_t2()`, `state()`, `pulse_shape()`, `kfigure()`, `dims()`, `volplot()`, `ktitle()`, `imaging()`, `scale_figure()`, `subplot()`, `mri_2d_plot()`, `apodisation()`, `fftshift()`.
