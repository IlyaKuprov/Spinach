# examples/imaging/gradient_image_1d.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/gradient_image_1d.m`
- Signature: `gradient_image_1d()`
- Total lines: 78

## Purpose

1D imaging experiment with a hard pulse in the presence of diffusion and flow. Calculation time: seconds. Ahmed Allami Ilya Kuprov

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Isotopes; implemented by `sys.isotopes={'1H'}`.
- Lines 14-15: Magnetic induction; implemented by `sys.magnet=5.9`.
- Lines 17-18: Chemical shifts; implemented by `inter.zeeman.scalar={1.0}`.
- Lines 20-21: Relaxation model; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 27-28: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 35-36: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 46-47: Sample geometry; implemented by `parameters.dims=0.30`.
- Lines 51-52: Relaxation phantom; implemented by `parameters.rlx_ph={zeros(parameters.npts,1)}`.
- Lines 55-56: Initial and detection state phantoms; implemented by `parameters.rho0_ph={ones(parameters.npts,1)}`.
- Lines 61-62: Diffusion and flow; implemented by `parameters.u=1e-2*ones(parameters.npts,1)`.
- Lines 65-66: Run the simulation; implemented by `fid=imaging(spin_system,@basic_1d_hard,parameters)`.
- Lines 68-69: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sqsin'}})`.
- Lines 71-72: Run the Fourier transform; implemented by `mri=real(fftshift(fft(ifftshift(fid))))`.
- Lines 74-75: Plotting; implemented by `kfigure(); plot_1d(spin_system,mri,parameters)`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H'}`.
- Lines 15: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.0}`.
- Lines 21: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 22: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 23: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 24: computes `inter.r1_rates` using `inter.r1_rates={0.5}`.
- Lines 25: computes `inter.r2_rates` using `inter.r2_rates={2.0}`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 36: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 37: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 38: computes `parameters.offset` using `parameters.offset=0`.
- Lines 39: computes `parameters.sweep` using `parameters.sweep=500000`.
- Lines 40: computes `parameters.npoints` using `parameters.npoints=128`.
- Lines 41: computes `parameters.zerofill` using `parameters.zerofill=128`.
- Lines 42: computes `parameters.axis_units` using `parameters.axis_units='kHz'`.

## Implementation structure

- 1D imaging experiment with a hard pulse in
- the presence of diffusion and flow.
- Calculation time: seconds.
- Ahmed Allami
- Ilya Kuprov
- Isotopes
- Magnetic induction
- Chemical shifts
- Relaxation model
- Basis set
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `state()`, `imaging()`, `apodisation()`, `fftshift()`, `ifftshift()`, `kfigure()`, `plot_1d()`.
