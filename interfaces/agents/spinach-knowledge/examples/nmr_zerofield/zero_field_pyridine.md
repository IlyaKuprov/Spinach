# examples/nmr_zerofield/zero_field_pyridine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_zerofield/zero_field_pyridine.m`
- Signature: `zero_field_pyridine()`
- Total lines: 66

## Purpose

Zero-field NMR spectroscopy -15N pyridine. Set to reproduce Figure 3 from http://dx.doi.org/10.1021/ja2112405 Calculation time: seconds

## Physical / mathematical content

- Zero- and ultralow-field NMR examples. The main physics is the crossover from Zeeman-dominated spectra to J-dominated spectra, with coherent evolution in near-zero field and detection of low-frequency transitions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnetic field; implemented by `sys.magnet=0`.
- Lines 13-14: Spin system; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','15N'}`.
- Lines 16-17: Interactions; implemented by `inter.coupling.scalar{1,2}= 4.88`.
- Lines 34-35: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 38-39: Sequence parameters; implemented by `parameters.sweep=60`.
- Lines 49-50: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 53-54: Simulation; implemented by `fid=liquid(spin_system,@zerofield,parameters,'labframe')`.
- Lines 56-57: Apodisation; implemented by `fid=apodisation(spin_system,fid-mean(fid),{{'exp',12}})`.
- Lines 59-60: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 62-63: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=0`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','15N'}`.
- Lines 17: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}= 4.88`.
- Lines 18: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}= 4.88`.
- Lines 19: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}= 0.97`.
- Lines 20: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}= 0.97`.
- Lines 21: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}= 1.83`.
- Lines 22: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}= 1.83`.
- Lines 23: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}= -0.12`.
- Lines 24: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}= 7.62`.
- Lines 25: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}= 7.62`.
- Lines 26: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}= 1.38`.
- Lines 27: computes `inter.coupling.scalar{1,6}` using `inter.coupling.scalar{1,6}= -10.93`.
- Lines 28: computes `inter.coupling.scalar{5,6}` using `inter.coupling.scalar{5,6}= -10.93`.
- Lines 29: computes `inter.coupling.scalar{2,6}` using `inter.coupling.scalar{2,6}= -1.47`.
- Lines 30: computes `inter.coupling.scalar{4,6}` using `inter.coupling.scalar{4,6}= -1.47`.
- Lines 31: computes `inter.coupling.scalar{3,6}` using `inter.coupling.scalar{3,6}= 0.27`.
- Lines 32: computes `inter.coupling.scalar{6,6}` using `inter.coupling.scalar{6,6}= 0.00`.

## Implementation structure

- Zero-field NMR spectroscopy -15N pyridine. Set to reproduce
- Figure 3 from http://dx.doi.org/10.1021/ja2112405
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
