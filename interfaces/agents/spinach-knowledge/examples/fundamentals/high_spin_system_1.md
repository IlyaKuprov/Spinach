# examples/fundamentals/high_spin_system_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/high_spin_system_1.m`
- Signature: `high_spin_system_1()`
- Total lines: 52

## Purpose

Pulse-acquire NMR spectrum in a system with a hypothetical scalar coupling to a 235U nucleus. The spectral lines should be split accordingly.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 12-13: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 16-17: Spin system; implemented by `sys.isotopes={'1H','235U','1H','1H'}`.
- Lines 23-24: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 27-28: Pulse sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 39-40: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 42-43: Apodization; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 45-46: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 48-49: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 14: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'1H','235U','1H','1H'}`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={-0.5 0.0 2.5 1.3}`.
- Lines 19: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=100`.
- Lines 20: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=50`.
- Lines 21: computes `inter.coupling.scalar{4,4}` using `inter.coupling.scalar{4,4}=0`.
- Lines 24: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 28: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 29: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 30: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 31: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 32: computes `parameters.offset` using `parameters.offset=0`.
- Lines 33: computes `parameters.sweep` using `parameters.sweep=3500`.
- Lines 34: computes `parameters.npoints` using `parameters.npoints=1024`.
- Lines 35: computes `parameters.zerofill` using `parameters.zerofill=4096`.
- Lines 36: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.

## Implementation structure

- Pulse-acquire NMR spectrum in a system with a hypothetical scalar
- coupling to a 235U nucleus. The spectral lines should be split
- accordingly.
- Magnet field
- Basis set
- Spin system
- Spinach housekeeping
- Pulse sequence parameters
- Simulation
- Apodization
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
