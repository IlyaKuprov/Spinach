# examples/nmr_liquids/pa_styrene.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/pa_styrene.m`
- Signature: `pa_styrene()`
- Total lines: 64

## Purpose

1H NMR spectrum of styrene, the calculation is performed in Hilbert space to demonstrate parallel propagation described in: Calculation time: seconds.

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Isotopes; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 15-16: Magnetic induction; implemented by `sys.magnet=7.046`.
- Lines 18-19: Chemical shifts; implemented by `inter.zeeman.scalar={6.720 6.720 6.587 6.587 6.527 6.063 5.085 4.555}`.
- Lines 21-22: Scalar couplings; implemented by `inter.coupling.scalar={ 0.0 1.9170 7.7884 0.6061 1.2447 -0.5347 0.0390 0.1582`.
- Lines 31-32: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 51-52: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 54-55: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'gauss',6}})`.
- Lines 57-58: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 60-61: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 16: computes `sys.magnet` using `sys.magnet=7.046`.
- Lines 19: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={6.720 6.720 6.587 6.587 6.527 6.063 5.085 4.555}`.
- Lines 22: computes `inter.coupling.scalar` using `inter.coupling.scalar={ 0.0 1.9170 7.7884 0.6061 1.2447 -0.5347 0.0390 0.1582`.
- Lines 32: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 33: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 41: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 42: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 43: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 44: computes `parameters.offset` using `parameters.offset=1700`.
- Lines 45: computes `parameters.sweep` using `parameters.sweep=1000`.
- Lines 46: computes `parameters.npoints` using `parameters.npoints=16384`.
- Lines 47: computes `parameters.zerofill` using `parameters.zerofill=65536`.
- Lines 48: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 49: computes `parameters.invert_axis` using `parameters.invert_axis=1`.
- Lines 52: computes `fid` using `fid=liquid(spin_system,@acquire,parameters,'nmr')`.

## Implementation structure

- 1H NMR spectrum of styrene, the calculation is performed in Hilbert
- space to demonstrate parallel propagation described in:
- Calculation time: seconds.
- Isotopes
- Magnetic induction
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
