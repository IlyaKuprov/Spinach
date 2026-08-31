# examples/extremes/phosphorus_cluster.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/extremes/phosphorus_cluster.m`
- Signature: `phosphorus_cluster()`
- Total lines: 112

## Purpose

Phosphorus system simulation for Gerhard Hagele. Done by brute force Liouville space time propagation. WARNING: needs 32+ CPU cores, 128+ GB of RAM and a strong FP64 capable Nvidia GPU. Run time on the above: hours

## Physical / mathematical content

- Extreme-regime examples. These scripts exercise Spinach in unusually large, stiff, high-field, low-field, or otherwise numerically demanding regimes where approximations, conditioning, and basis-size control are central.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Magnet induction; implemented by `sys.magnet=7.04`.
- Lines 16-20: Isotopes; implemented by `sys.isotopes={'31P','31P','31P','31P','31P','31P','31P', '1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H',…`.
- Lines 22-23: Chemical shifts; implemented by `inter.zeeman.scalar=cell(1,34)`.
- Lines 35-36: J-couplings; implemented by `inter.coupling.scalar=cell(34,34)`.
- Lines 68-69: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 76-77: Symmetry; implemented by `bas.sym_group={'S3','S3','S3'}`.
- Lines 80-81: Greedy parallelisation; implemented by `sys.enable={'greedy'}`.
- Lines 83-84: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 87-88: Sequence parameters; implemented by `parameters.spins={'31P'}`.
- Lines 99-100: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 102-103: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',5}})`.
- Lines 105-106: Fourier transform; implemented by `spectrum=real(fftshift(fft(fid,parameters.zerofill)))`.
- Lines 108-109: Plotting; implemented by `kfigure(); plot_1d(spin_system,spectrum,parameters)`.

### Control flow inferred from the code

- Line 31: `for` loop over `n=8:34`.
- Line 58: `for` loop over `n=8:16`.
- Line 61: `for` loop over `n=17:25`.
- Line 64: `for` loop over `n=26:34`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=7.04`.
- Lines 17-20: computes `sys.isotopes` using `sys.isotopes={'31P','31P','31P','31P','31P','31P','31P', '1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H',…`.
- Lines 23: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=cell(1,34)`.
- Lines 24: computes `inter.zeeman.scalar{1}` using `inter.zeeman.scalar{1}=-99.6`.
- Lines 25: computes `inter.zeeman.scalar{2}` using `inter.zeeman.scalar{2}=-0.33`.
- Lines 26: computes `inter.zeeman.scalar{3}` using `inter.zeeman.scalar{3}=-0.33`.
- Lines 27: computes `inter.zeeman.scalar{4}` using `inter.zeeman.scalar{4}=-0.33`.
- Lines 28: computes `inter.zeeman.scalar{5}` using `inter.zeeman.scalar{5}=-156.7`.
- Lines 29: computes `inter.zeeman.scalar{6}` using `inter.zeeman.scalar{6}=-156.7`.
- Lines 30: computes `inter.zeeman.scalar{7}` using `inter.zeeman.scalar{7}=-156.7`.
- Lines 32: computes `inter.zeeman.scalar{n}` using `inter.zeeman.scalar{n}=0.21`.
- Lines 36: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(34,34)`.
- Lines 37: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=-323.22`.
- Lines 38: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=-323.22`.
- Lines 39: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=-323.22`.
- Lines 40: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}= 46.18`.
- Lines 41: computes `inter.coupling.scalar{1,6}` using `inter.coupling.scalar{1,6}= 46.18`.
- Lines 42: computes `inter.coupling.scalar{1,7}` using `inter.coupling.scalar{1,7}= 46.18`.

## Implementation structure

- Phosphorus system simulation for Gerhard Hagele. Done by brute
- force Liouville space time propagation.
- WARNING: needs 32+ CPU cores, 128+ GB of RAM and
- a strong FP64 capable Nvidia GPU.
- Run time on the above: hours
- Magnet induction
- Isotopes
- Chemical shifts
- J-couplings
- Basis set
- Symmetry
- Greedy parallelisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
