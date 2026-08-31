# examples/extremes/high_symmetry_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/extremes/high_symmetry_2.m`
- Signature: `high_symmetry_2()`
- Total lines: 89

## Purpose

31P NMR spectrum of a large and highly symmetric spin system with two tert-butyl groups supplied by Eberhard Matern. Done by brute force time propagation in Hilbert space. WARNING: needs 32+ CPU cores and 128+ GB of RAM. Run time on the above: hours

## Physical / mathematical content

- Extreme-regime examples. These scripts exercise Spinach in unusually large, stiff, high-field, low-field, or otherwise numerically demanding regimes where approximations, conditioning, and basis-size control are central.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-17: Isotopes; implemented by `sys.isotopes={'31P','31P','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 19-20: Magnetic induction; implemented by `sys.magnet=9.39798`.
- Lines 22-28: Chemical shifts; implemented by `inter.zeeman.scalar={-43.844, -43.844, 4.090, 4.090, 1.354, 1.354, 1.354, 1.354, 1.354, 1.354 1.354, 1.354, 1.354, 1.354, 1.354, 1.354, 1.354, 1.354, 1.354 1.354, 1.354,…`.
- Lines 30-31: Scalar couplings; implemented by `inter.coupling.scalar=cell(22)`.
- Lines 56-57: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 60-61: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 64-65: Sequence parameters; implemented by `parameters.spins={'31P'}`.
- Lines 76-77: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 79-80: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',5}})`.
- Lines 82-83: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 85-86: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 15-17: computes `sys.isotopes` using `sys.isotopes={'31P','31P','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 20: computes `sys.magnet` using `sys.magnet=9.39798`.
- Lines 23-28: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={-43.844, -43.844, 4.090, 4.090, 1.354, 1.354, 1.354, 1.354, 1.354, 1.354 1.354, 1.354, 1.354, 1.354, 1.354, 1.354, 1.354, 1.354, 1.354 1.354, 1.354,…`.
- Lines 31: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(22)`.
- Lines 32: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}= 301.99`.
- Lines 33: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=-321.62`.
- Lines 34: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=-321.62`.
- Lines 35: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=-19.15`.
- Lines 36: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=-19.15`.
- Lines 37: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}= 15.63`.
- Lines 38: computes `inter.coupling.scalar{1,6}` using `inter.coupling.scalar{1,6}= 15.63`.
- Lines 39: computes `inter.coupling.scalar{1,7}` using `inter.coupling.scalar{1,7}= 15.63`.
- Lines 40: computes `inter.coupling.scalar{1,8}` using `inter.coupling.scalar{1,8}= 15.63`.
- Lines 41: computes `inter.coupling.scalar{1,9}` using `inter.coupling.scalar{1,9}= 15.63`.
- Lines 42: computes `inter.coupling.scalar{1,10}` using `inter.coupling.scalar{1,10}=15.63`.
- Lines 43: computes `inter.coupling.scalar{1,11}` using `inter.coupling.scalar{1,11}=15.63`.
- Lines 44: computes `inter.coupling.scalar{1,12}` using `inter.coupling.scalar{1,12}=15.63`.
- Lines 45: computes `inter.coupling.scalar{1,13}` using `inter.coupling.scalar{1,13}=15.63`.

## Implementation structure

- 31P NMR spectrum of a large and highly symmetric spin system
- with two tert-butyl groups supplied by Eberhard Matern. Done
- by brute force time propagation in Hilbert space.
- WARNING: needs 32+ CPU cores and 128+ GB of RAM.
- Run time on the above: hours
- Isotopes
- Magnetic induction
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
