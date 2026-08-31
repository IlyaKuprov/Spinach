# examples/nmr_zerofield/small_field_formic_acid.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_zerofield/small_field_formic_acid.m`
- Signature: `small_field_formic_acid()`
- Total lines: 52

## Purpose

Zero-field NMR spectroscopy -15N pyridine. Set to reproduce Fig 2 from http://dx.doi.org/10.1103/PhysRevLett.107.107601 Calculation time: seconds

## Physical / mathematical content

- Zero- and ultralow-field NMR examples. The main physics is the crossover from Zeeman-dominated spectra to J-dominated spectra, with coherent evolution in near-zero field and detection of low-frequency transitions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnetic field; implemented by `sys.magnet=1.76e-7`.
- Lines 13-14: Spin system; implemented by `sys.isotopes={'1H','13C'}`.
- Lines 16-17: Interactions; implemented by `inter.coupling.scalar=cell(2)`.
- Lines 20-21: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 24-25: Sequence parameters; implemented by `parameters.sweep=600`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Simulation; implemented by `fid=liquid(spin_system,@zerofield,parameters,'labframe')`.
- Lines 42-43: Apodisation; implemented by `fid=apodisation(spin_system,fid-mean(fid),{{'exp',12}})`.
- Lines 45-46: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 48-49: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=1.76e-7`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','13C'}`.
- Lines 17: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2)`.
- Lines 18: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=221`.
- Lines 21: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `parameters.sweep` using `parameters.sweep=600`.
- Lines 26: computes `parameters.npoints` using `parameters.npoints=8196`.
- Lines 27: computes `parameters.zerofill` using `parameters.zerofill=16384`.
- Lines 28: computes `parameters.offset` using `parameters.offset=0`.
- Lines 29: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 30: computes `parameters.axis_units` using `parameters.axis_units='Hz'`.
- Lines 31: computes `parameters.invert_axis` using `parameters.invert_axis=0`.
- Lines 32: computes `parameters.flip_angle` using `parameters.flip_angle=pi/2`.
- Lines 33: computes `parameters.detection` using `parameters.detection='uniaxial'`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `fid` using `fid=liquid(spin_system,@zerofield,parameters,'labframe')`.
- Lines 46: computes `spectrum` using `spectrum=fftshift(fft(fid,parameters.zerofill))`.

## Implementation structure

- Zero-field NMR spectroscopy -15N pyridine. Set to reproduce
- Fig 2 from http://dx.doi.org/10.1103/PhysRevLett.107.107601
- Calculation time: seconds
- Magnetic field
- Spin system
- Interactions
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- Fourier transform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
