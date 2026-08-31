# examples/fundamentals/strong_coupling.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/strong_coupling.m`
- Signature: `strong_coupling()`
- Total lines: 53

## Purpose

A garden variety strongly coupled two-spin system.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 7-8: Isotopes; implemented by `sys.isotopes={'1H','1H'}`.
- Lines 10-11: Magnetic induction; implemented by `sys.magnet=5.9`.
- Lines 13-14: Chemical shifts; implemented by `inter.zeeman.scalar={1.0 1.5}`.
- Lines 16-17: Scalar couplings; implemented by `inter.coupling.scalar=cell(2,2)`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 40-41: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 43-44: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',10}})`.
- Lines 46-47: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 49-50: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 8: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 11: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 14: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.0 1.5}`.
- Lines 17: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2,2)`.
- Lines 18: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=7.0`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 29: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 30: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 31: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 32: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 33: computes `parameters.offset` using `parameters.offset=300`.
- Lines 34: computes `parameters.sweep` using `parameters.sweep=300`.
- Lines 35: computes `parameters.npoints` using `parameters.npoints=1024`.
- Lines 36: computes `parameters.zerofill` using `parameters.zerofill=4096`.
- Lines 37: computes `parameters.axis_units` using `parameters.axis_units='Hz'`.
- Lines 38: computes `parameters.invert_axis` using `parameters.invert_axis=1`.

## Implementation structure

- A garden variety strongly coupled two-spin system.
- Isotopes
- Magnetic induction
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
