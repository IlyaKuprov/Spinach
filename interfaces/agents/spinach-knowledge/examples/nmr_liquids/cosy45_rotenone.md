# examples/nmr_liquids/cosy45_rotenone.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/cosy45_rotenone.m`
- Signature: `cosy45_rotenone()`
- Total lines: 83

## Purpose

Magnitude mode COSY-45 spectrum of rotenone. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-15: Spin system; implemented by `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H'}`.
- Lines 17-18: Interactions; implemented by `sys.magnet=5.9`.
- Lines 42-43: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 46-47: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 54-55: Sequence parameters; implemented by `parameters.angle=pi/4`.
- Lines 63-64: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 67-68: Simulation; implemented by `fid=liquid(spin_system,@cosy,parameters,'nmr')`.
- Lines 70-71: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'cos'},{'cos'}})`.
- Lines 73-75: Fourier transform; implemented by `spectrum=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 77-78: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 13-15: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H','1H','1H','1H','1H','1H', '1H','1H','1H','1H'}`.
- Lines 18: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 19-21: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={6.72 6.40 4.13 4.56 4.89 6.46 7.79 3.79 2.91 3.27 5.19 4.89 5.03 1.72 1.72 1.72 3.72 3.72 3.72 3.76 3.76 3.76}`.
- Lines 22: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=12.1`.
- Lines 23: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=3.1`.
- Lines 24: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=1.0`.
- Lines 25: computes `inter.coupling.scalar{3,8}` using `inter.coupling.scalar{3,8}=1.0`.
- Lines 26: computes `inter.coupling.scalar{1,8}` using `inter.coupling.scalar{1,8}=1.0`.
- Lines 27: computes `inter.coupling.scalar{6,7}` using `inter.coupling.scalar{6,7}=8.6`.
- Lines 28: computes `inter.coupling.scalar{5,8}` using `inter.coupling.scalar{5,8}=4.1`.
- Lines 29: computes `inter.coupling.scalar{7,9}` using `inter.coupling.scalar{7,9}=0.7`.
- Lines 30: computes `inter.coupling.scalar{7,10}` using `inter.coupling.scalar{7,10}=0.7`.
- Lines 31: computes `inter.coupling.scalar{9,10}` using `inter.coupling.scalar{9,10}=15.8`.
- Lines 32: computes `inter.coupling.scalar{10,11}` using `inter.coupling.scalar{10,11}=9.8`.
- Lines 33: computes `inter.coupling.scalar{9,11}` using `inter.coupling.scalar{9,11}=8.1`.
- Lines 34: computes `inter.coupling.scalar{13,14}` using `inter.coupling.scalar{13,14}=1.5`.
- Lines 35: computes `inter.coupling.scalar{12,14}` using `inter.coupling.scalar{12,14}=0.9`.
- Lines 36: computes `inter.coupling.scalar{13,15}` using `inter.coupling.scalar{13,15}=1.5`.

## Implementation structure

- Magnitude mode COSY-45 spectrum of rotenone.
- Calculation time: minutes
- Spin system
- Interactions
- Algorithmic options
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
