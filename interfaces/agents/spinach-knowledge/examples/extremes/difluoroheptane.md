# examples/extremes/difluoroheptane.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/extremes/difluoroheptane.m`
- Signature: `difluoroheptane()`
- Total lines: 122

## Purpose

19F NMR spectrum of anti-3,4-difluoroheptane (16 spins) by explicit time-domain evolution in Liouville space. WARNING: needs 32 CPU cores, 128 GB of RAM and a Titan V or later. Run time on the above: minutes

## Physical / mathematical content

- Extreme-regime examples. These scripts exercise Spinach in unusually large, stiff, high-field, low-field, or otherwise numerically demanding regimes where approximations, conditioning, and basis-size control are central.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnet induction; implemented by `sys.magnet=11.7464`.
- Lines 17-20: Isotopes; implemented by `sys.isotopes={'12C', '12C', '12C', '12C', '12C', '12C', '12C', '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H', '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H'}`.
- Lines 22-23: Shifts; implemented by `inter.zeeman.scalar=cell(1,23)`.
- Lines 41-42: J-couplings; implemented by `inter.coupling.scalar=cell(23)`.
- Lines 78-79: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 90-91: Greedy parallelisation; implemented by `sys.enable={'greedy'}`.
- Lines 93-94: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 97-98: Sequence parameters; implemented by `parameters.spins={'19F'}`.
- Lines 109-110: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 112-113: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 115-116: Fourier transform; implemented by `spectrum=real(fftshift(fft(fid,parameters.zerofill)))`.
- Lines 118-119: Plotting; implemented by `kfigure(); plot_1d(spin_system,spectrum,parameters)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=11.7464`.
- Lines 18-20: computes `sys.isotopes` using `sys.isotopes={'12C', '12C', '12C', '12C', '12C', '12C', '12C', '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H', '1H', '1H', '19F', '1H', '1H', '1H', '1H', '1H'}`.
- Lines 23: computes `inter.zeeman.scalar` using `inter.zeeman.scalar=cell(1,23)`.
- Lines 24: computes `inter.zeeman.scalar{14}` using `inter.zeeman.scalar{14}= 1.0092`.
- Lines 25: computes `inter.zeeman.scalar{15}` using `inter.zeeman.scalar{15}= 1.0092`.
- Lines 26: computes `inter.zeeman.scalar{16}` using `inter.zeeman.scalar{16}= 1.0092`.
- Lines 27: computes `inter.zeeman.scalar{21}` using `inter.zeeman.scalar{21}= 1.0092`.
- Lines 28: computes `inter.zeeman.scalar{22}` using `inter.zeeman.scalar{22}= 1.0092`.
- Lines 29: computes `inter.zeeman.scalar{23}` using `inter.zeeman.scalar{23}= 1.0092`.
- Lines 30: computes `inter.zeeman.scalar{11}` using `inter.zeeman.scalar{11}= 4.6834`.
- Lines 31: computes `inter.zeeman.scalar{17}` using `inter.zeeman.scalar{17}= 4.6834`.
- Lines 32: computes `inter.zeeman.scalar{10}` using `inter.zeeman.scalar{10}=-184.1865`.
- Lines 33: computes `inter.zeeman.scalar{18}` using `inter.zeeman.scalar{18}=-184.1865`.
- Lines 34: computes `inter.zeeman.scalar{8}` using `inter.zeeman.scalar{8}= 1.7970`.
- Lines 35: computes `inter.zeeman.scalar{9}` using `inter.zeeman.scalar{9}= 1.7970`.
- Lines 36: computes `inter.zeeman.scalar{13}` using `inter.zeeman.scalar{13}= 1.6942`.
- Lines 37: computes `inter.zeeman.scalar{20}` using `inter.zeeman.scalar{20}= 1.6942`.
- Lines 38: computes `inter.zeeman.scalar{19}` using `inter.zeeman.scalar{19}= 1.6370`.

## Implementation structure

- 19F NMR spectrum of anti-3,4-difluoroheptane (16 spins) by
- explicit time-domain evolution in Liouville space.
- WARNING: needs 32 CPU cores, 128 GB of RAM and
- a Titan V or later.
- Run time on the above: minutes
- Magnet induction
- Isotopes
- Shifts
- J-couplings
- Basis set
- Greedy parallelisation
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `false()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
