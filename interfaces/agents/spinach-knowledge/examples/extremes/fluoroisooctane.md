# examples/extremes/fluoroisooctane.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/extremes/fluoroisooctane.m`
- Signature: `fluoroisooctane()`
- Total lines: 96

## Purpose

A deliberately adversarial example from Art Bochevarov at Schodinger Inc. In this case, IK-2 approximation in Liou- ville space generates an exceedingly large basis set; the calculation must instead be performed in Hilbert space with permutation symmetry factorisation. Calculation time: hours.

## Physical / mathematical content

- Extreme-regime examples. These scripts exercise Spinach in unusually large, stiff, high-field, low-field, or otherwise numerically demanding regimes where approximations, conditioning, and basis-size control are central.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnet induction; implemented by `sys.magnet=11.74`.
- Lines 17-20: Isotopes; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','19F','1H', '1H','1H','1H','1H','1H','1H'}`.
- Lines 22-28: Chemical shifts; implemented by `inter.zeeman.scalar={0.994, 0.994, 0.994, 0.994, 0.994, 0.994, 0.994, 0.994, 0.994, 4.165, 0.000, 1.938, 1.039, 1.039, 1.039, 1.062, 1.062, 1.062}`.
- Lines 30-31: Larger J-couplings; implemented by `inter.coupling.scalar=cell(18,18)`.
- Lines 42-43: Smaller J-couplings, tert-butyl; implemented by `inter.coupling.scalar{1,11}=1.0`.
- Lines 53-54: Smaller J-couplings, isopropyl; implemented by `inter.coupling.scalar{13,11}=1.0`.
- Lines 61-62: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 67-68: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 71-72: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 83-84: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 86-87: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'gauss',6}})`.
- Lines 89-90: Fourier transform; implemented by `spectrum=real(fftshift(fft(fid,parameters.zerofill)))`.
- Lines 92-93: Plotting; implemented by `plot_1d(spin_system,spectrum,parameters)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=11.74`.
- Lines 18-20: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','19F','1H', '1H','1H','1H','1H','1H','1H'}`.
- Lines 23-28: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.994, 0.994, 0.994, 0.994, 0.994, 0.994, 0.994, 0.994, 0.994, 4.165, 0.000, 1.938, 1.039, 1.039, 1.039, 1.062, 1.062, 1.062}`.
- Lines 31: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(18,18)`.
- Lines 32: computes `inter.coupling.scalar{10,11}` using `inter.coupling.scalar{10,11}=48.2`.
- Lines 33: computes `inter.coupling.scalar{11,12}` using `inter.coupling.scalar{11,12}=23.6`.
- Lines 34: computes `inter.coupling.scalar{10,12}` using `inter.coupling.scalar{10,12}= 6.5`.
- Lines 35: computes `inter.coupling.scalar{12,13}` using `inter.coupling.scalar{12,13}= 6.0`.
- Lines 36: computes `inter.coupling.scalar{12,14}` using `inter.coupling.scalar{12,14}= 6.0`.
- Lines 37: computes `inter.coupling.scalar{12,15}` using `inter.coupling.scalar{12,15}= 6.0`.
- Lines 38: computes `inter.coupling.scalar{12,16}` using `inter.coupling.scalar{12,16}= 6.0`.
- Lines 39: computes `inter.coupling.scalar{12,17}` using `inter.coupling.scalar{12,17}= 6.0`.
- Lines 40: computes `inter.coupling.scalar{12,18}` using `inter.coupling.scalar{12,18}= 6.0`.
- Lines 43: computes `inter.coupling.scalar{1,11}` using `inter.coupling.scalar{1,11}=1.0`.
- Lines 44: computes `inter.coupling.scalar{2,11}` using `inter.coupling.scalar{2,11}=1.0`.
- Lines 45: computes `inter.coupling.scalar{3,11}` using `inter.coupling.scalar{3,11}=1.0`.
- Lines 46: computes `inter.coupling.scalar{4,11}` using `inter.coupling.scalar{4,11}=1.0`.
- Lines 47: computes `inter.coupling.scalar{5,11}` using `inter.coupling.scalar{5,11}=1.0`.

## Implementation structure

- A deliberately adversarial example from Art Bochevarov at
- Schodinger Inc. In this case, IK-2 approximation in Liou-
- ville space generates an exceedingly large basis set; the
- calculation must instead be performed in Hilbert space
- with permutation symmetry factorisation.
- Calculation time: hours.
- Magnet induction
- Isotopes
- Chemical shifts
- Larger J-couplings
- Smaller J-couplings, tert-butyl
- Smaller J-couplings, isopropyl

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `plot_1d()`.
