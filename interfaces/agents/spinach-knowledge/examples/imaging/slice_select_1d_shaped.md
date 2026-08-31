# examples/imaging/slice_select_1d_shaped.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/slice_select_1d_shaped.m`
- Signature: `slice_select_1d_shaped()`
- Total lines: 88

## Purpose

Slice selection example using a one-dimensional sample and a shaped slice selection pulse in the presence of diffusion and flow. Calculation time: seconds. Ahmed Allami Ilya Kuprov

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Isotopes; implemented by `sys.isotopes={'1H'}`.
- Lines 15-16: Magnetic induction; implemented by `sys.magnet=5.9`.
- Lines 18-19: Chemical shifts; implemented by `inter.zeeman.scalar={0.0}`.
- Lines 21-22: Relaxation model; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 28-29: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 47-48: Pulse parameters; implemented by `pulse_nsteps=50; pulse_time=0.5e-4`.
- Lines 56-57: Sample geometry; implemented by `parameters.dims=0.30`.
- Lines 61-62: Relaxation phantom; implemented by `parameters.rlx_ph={zeros(parameters.npts,1)}`.
- Lines 65-66: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 71-72: Diffusion and flow; implemented by `parameters.u=1e-2*ones(parameters.npts,1)`.
- Lines 75-76: Run the pulse sequence in the imaging context; implemented by `fid=imaging(spin_system,@slice_select_1d,parameters)`.
- Lines 78-79: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sqsin'}})`.
- Lines 81-82: Run the Fourier transform; implemented by `mri=real(fftshift(fft(ifftshift(fid))))`.
- Lines 84-85: Plotting; implemented by `kfigure(); plot_1d(spin_system,mri,parameters)`.

### Key state/data transformations

- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 16: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 19: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0}`.
- Lines 22: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 23: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 24: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 25: computes `inter.r1_rates` using `inter.r1_rates={30}`.
- Lines 26: computes `inter.r2_rates` using `inter.r2_rates={70}`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 38: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 39: computes `parameters.offset` using `parameters.offset=0`.
- Lines 40: computes `parameters.sweep` using `parameters.sweep=500000`.
- Lines 41: computes `parameters.npoints` using `parameters.npoints=128`.
- Lines 42: computes `parameters.axis_units` using `parameters.axis_units='kHz'`.
- Lines 43: computes `parameters.invert_axis` using `parameters.invert_axis=1`.

## Implementation structure

- Slice selection example using a one-dimensional sample and
- a shaped slice selection pulse in the presence of diffusion
- and flow.
- Calculation time: seconds.
- Ahmed Allami
- Ilya Kuprov
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation model
- Basis set
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `pulse_shape()`, `relaxation()`, `state()`, `imaging()`, `apodisation()`, `fftshift()`, `ifftshift()`, `kfigure()`, `plot_1d()`.
