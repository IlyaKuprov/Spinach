# examples/nmr_zerofield/zero_field_methanol.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_zerofield/zero_field_methanol.m`
- Signature: `zero_field_methanol()`
- Total lines: 55

## Purpose

Zero-field NMR spectroscopy -13C methanol. Set to reproduce Figure 1 from http://dx.doi.org/10.1016/j.cplett.2013.06.042 Calculation time: seconds

## Physical / mathematical content

- Zero- and ultralow-field NMR examples. The main physics is the crossover from Zeeman-dominated spectra to J-dominated spectra, with coherent evolution in near-zero field and detection of low-frequency transitions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Spin system; implemented by `sys.isotopes={'1H','1H','1H','13C'}`.
- Lines 13-14: Interactions; implemented by `sys.magnet=0`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-28: Sequence parameters; implemented by `parameters.sweep=700`.
- Lines 38-39: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Simulation; implemented by `fid=liquid(spin_system,@zerofield,parameters,'labframe')`.
- Lines 45-46: Apodisation; implemented by `fid=apodisation(spin_system,fid-mean(fid),{{'exp',6}})`.
- Lines 48-49: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 51-52: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','13C'}`.
- Lines 14: computes `sys.magnet` using `sys.magnet=0`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0 0 0 0}`.
- Lines 16: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=141.0`.
- Lines 17: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=141.0`.
- Lines 18: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=141.0`.
- Lines 19: computes `inter.coupling.scalar{4,4}` using `inter.coupling.scalar{4,4}=0`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `bas.sym_group` using `bas.sym_group={'S3'}`.
- Lines 25: computes `bas.sym_spins` using `bas.sym_spins={[1 2 3]}`.
- Lines 28: computes `parameters.sweep` using `parameters.sweep=700`.
- Lines 29: computes `parameters.npoints` using `parameters.npoints=4096`.
- Lines 30: computes `parameters.zerofill` using `parameters.zerofill=16384`.
- Lines 31: computes `parameters.offset` using `parameters.offset=0`.
- Lines 32: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 33: computes `parameters.axis_units` using `parameters.axis_units='Hz'`.
- Lines 34: computes `parameters.invert_axis` using `parameters.invert_axis=0`.

## Implementation structure

- Zero-field NMR spectroscopy -13C methanol. Set to reproduce
- Figure 1 from http://dx.doi.org/10.1016/j.cplett.2013.06.042
- Calculation time: seconds
- Spin system
- Interactions
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
