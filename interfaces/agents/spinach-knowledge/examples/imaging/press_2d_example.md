# examples/imaging/press_2d_example.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/imaging/press_2d_example.m`
- Signature: `press_2d_example()`
- Total lines: 113

## Purpose

2D PRESS example. Three independent spin systems are localised in three spots of a 2D sample. The spots are slectively excited and their NMR spectra recorded. The followng are the frequencies to excite the three substances: parameters.rf_frq_list={-120e3 -100e3} -substance A parameters.rf_frq_list={-80e3 -10e3} -substance B parameters.rf_frq_list={+30e3 +100e3} -substance C Simulation time: minutes, faster with a Tes

## Physical / mathematical content

- MRI and spectroscopic-imaging examples. These files combine gradient terms, spatial encoding, diffusion, slice selection, k-space sampling, and Fourier reconstruction, generally within Fokker-Planck or explicit spatial-grid descriptions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Magnetic induction; implemented by `sys.magnet=3.0`.
- Lines 22-23: Spin systems; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Disable path tracing; implemented by `sys.disable={'pt'}`.
- Lines 39-43: This needs a GPU sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 46-47: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 59-60: Pulse parameters; implemented by `parameters.rf_phi={pi/2 pi/2}`.
- Lines 66-67: Sample geometry; implemented by `parameters.dims=[0.30 0.25]`.
- Lines 71-72: Relaxation phantom; implemented by `parameters.rlx_ph={zeros(parameters.npts)}`.
- Lines 75-76: Initial and detection state phantoms; implemented by `[X,Y]=ndgrid(1:parameters.npts(1),1:parameters.npts(2))`.
- Lines 86-87: Show the phantom; implemented by `kfigure(); scale_figure([2.0 1.0])`.
- Lines 94-95: Run active volume diagnostics; implemented by `phan=imaging(spin_system,@press_voxel_2d,parameters)`.
- Lines 99-100: Run PRESS 2D; implemented by `fid=imaging(spin_system,@press_2d,parameters)`.
- Lines 102-103: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sqcos'}})`.
- Lines 105-106: Magnitude Fourier transform; implemented by `spec=abs(fftshift(fft(ifftshift(fid))))`.
- Lines 108-109: Plotting; implemented by `subplot(1,3,3); plot_1d(spin_system,spec,parameters)`.

### Key state/data transformations

- Lines 20: computes `sys.magnet` using `sys.magnet=3.0`.
- Lines 23: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H'}`.
- Lines 24: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={-3,3,-2,2,-1,1}`.
- Lines 25: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(6)`.
- Lines 26: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=10`.
- Lines 27: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=20`.
- Lines 28: computes `inter.coupling.scalar{5,6}` using `inter.coupling.scalar{5,6}=30`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 33: computes `bas.space_level` using `bas.space_level=1`.
- Lines 34: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 37: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 43: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 47: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 48: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 49: computes `parameters.offset` using `parameters.offset=0`.
- Lines 50: computes `parameters.sweep` using `parameters.sweep=1000`.
- Lines 51: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.

## Implementation structure

- 2D PRESS example. Three independent spin systems are localised
- in three spots of a 2D sample. The spots are slectively excited
- and their NMR spectra recorded. The followng are the frequencies
- to excite the three substances:
- parameters.rf_frq_list={-120e3 -100e3} -substance A
- parameters.rf_frq_list={-80e3 -10e3} -substance B
- parameters.rf_frq_list={+30e3 +100e3} -substance C
- Simulation time: minutes, faster with a Tesla V100 GPU.
- Magnetic induction
- Spin systems
- Basis set
- Disable path tracing

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `relaxation()`, `state()`, `kfigure()`, `scale_figure()`, `subplot()`, `mri_2d_plot()`, `ktitle()`, `imaging()`, `apodisation()`, `fftshift()`, `ifftshift()`, `plot_1d()`.
