# examples/fundamentals/symmetry_3.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/symmetry_3.m`
- Signature: `symmetry_3()`
- Total lines: 57

## Purpose

1H NMR spectrum of valine. Uses the fully symmetric irreducible representation of the S3(x)S3 group.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: System specificaton; implemented by `sys.magnet=11.7`.
- Lines 22-23: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 44-45: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 47-48: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',5}})`.
- Lines 50-51: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 53-54: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 9: computes `sys.magnet` using `sys.magnet=11.7`.
- Lines 10: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 11-12: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={3.5950 2.2580 1.0270 1.0270 1.0270 0.9760 0.9760 0.9760}`.
- Lines 13: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=4.34`.
- Lines 14: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=7.00`.
- Lines 15: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=7.00`.
- Lines 16: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}=7.00`.
- Lines 17: computes `inter.coupling.scalar{2,6}` using `inter.coupling.scalar{2,6}=7.00`.
- Lines 18: computes `inter.coupling.scalar{2,7}` using `inter.coupling.scalar{2,7}=7.00`.
- Lines 19: computes `inter.coupling.scalar{2,8}` using `inter.coupling.scalar{2,8}=7.00`.
- Lines 20: computes `inter.coupling.scalar{8,8}` using `inter.coupling.scalar{8,8}=0.00`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `bas.sym_spins` using `bas.sym_spins={[3 4 5],[6 7 8]}`.
- Lines 26: computes `bas.sym_group` using `bas.sym_group={'S3','S3'}`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 34: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.

## Implementation structure

- 1H NMR spectrum of valine. Uses the fully symmetric
- irreducible representation of the S3(x)S3 group.
- System specificaton
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
