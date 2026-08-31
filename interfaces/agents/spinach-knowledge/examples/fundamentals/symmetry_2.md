# examples/fundamentals/symmetry_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/symmetry_2.m`
- Signature: `symmetry_2()`
- Total lines: 68

## Purpose

Pulse-acquire NMR spectrum of a highly symmetric spin system provided by Andres Castillo. Uses the fully sym- metric irreducible representation of S3(x)S3(x)S3 per- mutation symmetry group.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Spin system specification; implemented by `sys.magnet=9.4`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 39-40: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 43-44: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 55-56: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 58-59: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 61-62: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 64-65: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 14: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.89 0.89 0.89 0.895 0.895 0.895 1.16 1.16 1.16 1.2 1.39 1.71 3.85}`.
- Lines 15-16: computes `inter.coupling.scalar` using `inter.coupling.scalar=num2cell( [ 0 0 0 0 0 0 0 0 0 0 0 6.7200 0`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 33: computes `bas.sym_group` using `bas.sym_group={'S3','S3','S3'}`.
- Lines 34: computes `bas.sym_spins` using `bas.sym_spins={[1 2 3],[4 5 6],[7 8 9]}`.
- Lines 35: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 36: computes `bas.space_level` using `bas.space_level=1`.
- Lines 37: computes `bas.projections` using `bas.projections=+1`.
- Lines 40: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 44: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 45: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 46: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 47: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 48: computes `parameters.offset` using `parameters.offset=800`.
- Lines 49: computes `parameters.sweep` using `parameters.sweep=2000`.

## Implementation structure

- Pulse-acquire NMR spectrum of a highly symmetric spin
- system provided by Andres Castillo. Uses the fully sym-
- metric irreducible representation of S3(x)S3(x)S3 per-
- mutation symmetry group.
- Spin system specification
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
