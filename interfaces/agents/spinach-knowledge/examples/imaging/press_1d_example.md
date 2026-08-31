# examples/imaging/press_1d_example.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/press_1d_example.m`
- Signature: `press_1d_example()`
- Total lines: 115

## Purpose

1D PRESS example. Three independent spin systems are localised in three areas of a 1D sample. The areas are slectively excited and their NMR spectra recorded. The followng are the frequencies to excite the three substances: pulse_frq=+100e3 -substances B and C pulse_frq=0; -substance C pulse_frq=-100e3 -substances A and C Simulation time: seconds, faster with a Tesla V100 GPU.

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Magnetic induction; implemented by `sys.magnet=3.0`.
- Lines 28-29: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Disable path tracing; implemented by `sys.disable={'pt'}`.
- Lines 37-41: This needs a GPU sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 55-56: Pulse parameters; implemented by `parameters.rf_phi=pi/2`.
- Lines 62-63: Sample geometry; implemented by `parameters.dims=0.30`.
- Lines 67-68: Relaxation phantom; implemented by `parameters.rlx_ph={zeros(parameters.npts,1)}`.
- Lines 71-72: Initial and detection state phantoms; implemented by `parameters.rho0_ph={[zeros(10,1); ones(30,1); zeros(60,1)]`.
- Lines 81-82: Show the phantom; implemented by `kfigure(); scale_figure([1.50 0.75])`.
- Lines 90-91: Run voxel selection diagnostics; implemented by `phan=imaging(spin_system,@press_voxel_1d,parameters)`.
- Lines 93-94: Plotting; implemented by `subplot(1,3,2)`.
- Lines 99-100: Run PRESS 1D; implemented by `parameters.sweep=1000`.
- Lines 104-105: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sqcos'}})`.
- Lines 107-108: Magnitude Fourier transform; implemented by `spec=abs(fftshift(fft(ifftshift(fid))))`.
- Lines 110-111: Plotting; implemented by `subplot(1,3,3); plot_1d(spin_system,spec,parameters)`.

### Key state/data transformations

- Lines 20: computes `sys.magnet` using `sys.magnet=3.0`.
- Lines 21: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 22: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={-3,3,-2,2,-1,1}`.
- Lines 23: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(6)`.
- Lines 24: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=10`.
- Lines 25: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=20`.
- Lines 26: computes `inter.coupling.scalar{5,6}` using `inter.coupling.scalar{5,6}=30`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 31: computes `bas.space_level` using `bas.space_level=1`.
- Lines 32: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 35: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 41: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 45: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 46: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 47: computes `parameters.offset` using `parameters.offset=0`.
- Lines 48: computes `parameters.sweep` using `parameters.sweep=500000`.
- Lines 49: computes `parameters.npoints` using `parameters.npoints=128`.

## Implementation structure

- 1D PRESS example. Three independent spin systems are localised
- in three areas of a 1D sample. The areas are slectively excited
- and their NMR spectra recorded. The followng are the frequencies
- to excite the three substances:
- pulse_frq=+100e3 -substances B and C
- pulse_frq=0; -substance C
- pulse_frq=-100e3 -substances A and C
- Simulation time: seconds, faster with a Tesla V100 GPU.
- Magnetic induction
- Basis set
- Disable path tracing
- This needs a GPU

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `state()`, `kfigure()`, `scale_figure()`, `subplot()`, `ylim()`, `ktitle()`, `kxlabel()`, `imaging()`, `apodisation()`, `fftshift()`, `ifftshift()`, `plot_1d()`.
