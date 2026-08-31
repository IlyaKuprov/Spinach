# examples/nmr_liquids/pa_naphtopyranone.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/pa_naphtopyranone.m`
- Signature: `pa_naphtopyranone()`
- Total lines: 73

## Purpose

NMR spectrum of 3-phenylmethylene-1H,3H-naphtho-[1,8-c,d]-pyran-1-one, magnetic parameters from: Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Magnetic induction; implemented by `sys.magnet=14.095`.
- Lines 16-17: Spin system; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 19-21: Chemical shifts; implemented by `inter.zeeman.scalar={8.345,7.741,8.097,8.354,7.784,8.330,7.059, 7.941,7.466,7.326,7.466,7.941}`.
- Lines 23-24: Scalar couplings; implemented by `inter.coupling.scalar=cell(12,12)`.
- Lines 38-39: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 44-45: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 48-49: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 60-61: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 63-64: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',10}})`.
- Lines 66-67: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 69-70: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=14.095`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H','1H'}`.
- Lines 20-21: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={8.345,7.741,8.097,8.354,7.784,8.330,7.059, 7.941,7.466,7.326,7.466,7.941}`.
- Lines 24: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(12,12)`.
- Lines 25: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=7.8`.
- Lines 26: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=0.9`.
- Lines 27: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=7.8`.
- Lines 28: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=8.4`.
- Lines 29: computes `inter.coupling.scalar{4,6}` using `inter.coupling.scalar{4,6}=1.2`.
- Lines 30: computes `inter.coupling.scalar{5,6}` using `inter.coupling.scalar{5,6}=7.2`.
- Lines 31: computes `inter.coupling.scalar{8,9}` using `inter.coupling.scalar{8,9}=7.8`.
- Lines 32: computes `inter.coupling.scalar{8,10}` using `inter.coupling.scalar{8,10}=1.2`.
- Lines 33: computes `inter.coupling.scalar{9,10}` using `inter.coupling.scalar{9,10}=7.8`.
- Lines 34: computes `inter.coupling.scalar{10,11}` using `inter.coupling.scalar{10,11}=7.8`.
- Lines 35: computes `inter.coupling.scalar{10,12}` using `inter.coupling.scalar{10,12}=1.2`.
- Lines 36: computes `inter.coupling.scalar{11,12}` using `inter.coupling.scalar{11,12}=7.8`.
- Lines 39: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 40: computes `bas.approximation` using `bas.approximation='IK-2'`.

## Implementation structure

- NMR spectrum of 3-phenylmethylene-1H,3H-naphtho-[1,8-c,d]-pyran-1-one,
- magnetic parameters from:
- Calculation time: seconds
- Magnetic induction
- Spin system
- Chemical shifts
- Scalar couplings
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
