# examples/nmr_zerofield/small_field_acetonitrile.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_zerofield/small_field_acetonitrile.m`
- Signature: `small_field_acetonitrile()`
- Total lines: 59

## Purpose

Small-field NMR spectroscopy -acetonitrile with 13C on the methyl group. Set to reproduce Figure 3 from Calculation time: seconds

## Physical / mathematical content

- Zero- and ultralow-field NMR examples. The main physics is the crossover from Zeeman-dominated spectra to J-dominated spectra, with coherent evolution in near-zero field and detection of low-frequency transitions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnetic field, 2.64 mG; implemented by `sys.magnet=2.64e-3*1e-4`.
- Lines 15-16: Spin system; implemented by `sys.isotopes={'1H','1H','1H','13C'}`.
- Lines 18-19: Interactions; implemented by `inter.coupling.scalar=cell(4,4)`.
- Lines 24-25: Temperature; implemented by `inter.temperature=298`.
- Lines 27-28: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 31-32: Sequence parameters; implemented by `parameters.sweep=700`.
- Lines 42-43: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 46-47: Simulation; implemented by `fid=liquid(spin_system,@zerofield,parameters,'labframe')`.
- Lines 49-50: Apodisation; implemented by `fid=apodisation(spin_system,fid-mean(fid),{{'exp',6}})`.
- Lines 52-53: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 55-56: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=2.64e-3*1e-4`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','13C'}`.
- Lines 19: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(4,4)`.
- Lines 20: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=136.200`.
- Lines 21: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=136.200`.
- Lines 22: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=136.200`.
- Lines 25: computes `inter.temperature` using `inter.temperature=298`.
- Lines 28: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `parameters.sweep` using `parameters.sweep=700`.
- Lines 33: computes `parameters.npoints` using `parameters.npoints=4096`.
- Lines 34: computes `parameters.zerofill` using `parameters.zerofill=16384`.
- Lines 35: computes `parameters.offset` using `parameters.offset=0`.
- Lines 36: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 37: computes `parameters.axis_units` using `parameters.axis_units='Hz'`.
- Lines 38: computes `parameters.invert_axis` using `parameters.invert_axis=0`.
- Lines 39: computes `parameters.flip_angle` using `parameters.flip_angle=pi/2`.
- Lines 40: computes `parameters.detection` using `parameters.detection='uniaxial'`.

## Implementation structure

- Small-field NMR spectroscopy -acetonitrile with 13C on the
- methyl group. Set to reproduce Figure 3 from
- Calculation time: seconds
- Magnetic field, 2.64 mG
- Spin system
- Interactions
- Temperature
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
