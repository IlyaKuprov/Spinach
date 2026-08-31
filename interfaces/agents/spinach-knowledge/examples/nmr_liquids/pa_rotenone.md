# examples/nmr_liquids/pa_rotenone.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/pa_rotenone.m`
- Signature: `pa_rotenone()`
- Total lines: 91

## Purpose

1H NMR spectrum of rotenone using T1/T2 relaxation model, magnetic parameters from: Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-16: Isotopes; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H'}`.
- Lines 18-19: Magnetic induction; implemented by `sys.magnet=5.9`.
- Lines 21-24: Chemical shifts; implemented by `inter.zeeman.scalar={6.72 6.40 4.13 4.56 4.89 6.46 7.79 3.79 2.91 3.27 5.19 4.89 5.03 1.72 1.72 1.72 3.72 3.72 3.72 3.76 3.76 3.76}`.
- Lines 26-27: Scalar couplings; implemented by `inter.coupling.scalar{3,4}=12.1`.
- Lines 47-48: Relaxation model; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 54-55: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 62-63: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 66-67: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 78-79: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 81-82: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'gauss',10}})`.
- Lines 84-85: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 87-88: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 14-16: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H'}`.
- Lines 19: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 22-24: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={6.72 6.40 4.13 4.56 4.89 6.46 7.79 3.79 2.91 3.27 5.19 4.89 5.03 1.72 1.72 1.72 3.72 3.72 3.72 3.76 3.76 3.76}`.
- Lines 27: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=12.1`.
- Lines 28: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=3.1`.
- Lines 29: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=1.0`.
- Lines 30: computes `inter.coupling.scalar{3,8}` using `inter.coupling.scalar{3,8}=1.0`.
- Lines 31: computes `inter.coupling.scalar{1,8}` using `inter.coupling.scalar{1,8}=1.0`.
- Lines 32: computes `inter.coupling.scalar{6,7}` using `inter.coupling.scalar{6,7}=8.6`.
- Lines 33: computes `inter.coupling.scalar{5,8}` using `inter.coupling.scalar{5,8}=4.1`.
- Lines 34: computes `inter.coupling.scalar{7,9}` using `inter.coupling.scalar{7,9}=0.7`.
- Lines 35: computes `inter.coupling.scalar{7,10}` using `inter.coupling.scalar{7,10}=0.7`.
- Lines 36: computes `inter.coupling.scalar{9,10}` using `inter.coupling.scalar{9,10}=15.8`.
- Lines 37: computes `inter.coupling.scalar{10,11}` using `inter.coupling.scalar{10,11}=9.8`.
- Lines 38: computes `inter.coupling.scalar{9,11}` using `inter.coupling.scalar{9,11}=8.1`.
- Lines 39: computes `inter.coupling.scalar{13,14}` using `inter.coupling.scalar{13,14}=1.5`.
- Lines 40: computes `inter.coupling.scalar{12,14}` using `inter.coupling.scalar{12,14}=0.9`.
- Lines 41: computes `inter.coupling.scalar{13,15}` using `inter.coupling.scalar{13,15}=1.5`.

## Implementation structure

- 1H NMR spectrum of rotenone using T1/T2 relaxation model,
- magnetic parameters from:
- Calculation time: seconds
- Isotopes
- Magnetic induction
- Chemical shifts
- Scalar couplings
- Relaxation model
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `num2cell()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
