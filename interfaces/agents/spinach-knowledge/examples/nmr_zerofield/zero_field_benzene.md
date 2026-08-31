# examples/nmr_zerofield/zero_field_benzene.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_zerofield/zero_field_benzene.m`
- Signature: `zero_field_benzene()`
- Total lines: 74

## Purpose

Zero-field NMR spectroscopy -benzene with one 13C nucleus. Set to reproduce Figure 2 from Calculation time: seconds

## Physical / mathematical content

- Zero- and ultralow-field NMR examples. The main physics is the crossover from Zeeman-dominated spectra to J-dominated spectra, with coherent evolution in near-zero field and detection of low-frequency transitions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnetic field; implemented by `sys.magnet=0`.
- Lines 15-16: Spin system; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H','13C'}`.
- Lines 18-19: Interactions; implemented by `inter.coupling.scalar{1,7}=158.363`.
- Lines 42-43: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 46-47: Sequence parameters; implemented by `parameters.sweep=400`.
- Lines 57-58: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 61-62: Simulation; implemented by `fid=liquid(spin_system,@zerofield,parameters,'labframe')`.
- Lines 64-65: Apodisation; implemented by `fid=apodisation(spin_system,fid-mean(fid),{{'exp',6}})`.
- Lines 67-68: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 70-71: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=0`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','13C'}`.
- Lines 19: computes `inter.coupling.scalar{1,7}` using `inter.coupling.scalar{1,7}=158.363`.
- Lines 20: computes `inter.coupling.scalar{2,7}` using `inter.coupling.scalar{2,7}=1.136`.
- Lines 21: computes `inter.coupling.scalar{3,7}` using `inter.coupling.scalar{3,7}=7.609`.
- Lines 22: computes `inter.coupling.scalar{4,7}` using `inter.coupling.scalar{4,7}=-1.285`.
- Lines 23: computes `inter.coupling.scalar{5,7}` using `inter.coupling.scalar{5,7}=7.609`.
- Lines 24: computes `inter.coupling.scalar{6,7}` using `inter.coupling.scalar{6,7}=1.136`.
- Lines 25: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=7.534`.
- Lines 26: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=1.381`.
- Lines 27: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=0.658`.
- Lines 28: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}=1.381`.
- Lines 29: computes `inter.coupling.scalar{1,6}` using `inter.coupling.scalar{1,6}=7.534`.
- Lines 30: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=7.543`.
- Lines 31: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=1.382`.
- Lines 32: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}=0.660`.
- Lines 33: computes `inter.coupling.scalar{2,6}` using `inter.coupling.scalar{2,6}=1.384`.
- Lines 34: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=7.543`.

## Implementation structure

- Zero-field NMR spectroscopy -benzene with one 13C
- nucleus. Set to reproduce Figure 2 from
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
