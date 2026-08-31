# examples/nmr_liquids/coloc_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/coloc_test.m`
- Signature: `coloc_test()`
- Total lines: 54

## Purpose

A simple COLOC pulse sequence example for a two-spin 1H-13C system with a long-range J-coupling. Calculation time: seconds.

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnet field; implemented by `sys.magnet=11.7`.
- Lines 14-15: Spin system; implemented by `sys.isotopes={'1H','13C'}`.
- Lines 17-18: Interactions; implemented by `inter.zeeman.scalar={4.0 75.0}`.
- Lines 22-23: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 26-27: Sequence parameters; implemented by `parameters.spins={'1H','13C'}`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Simulation; implemented by `fid=liquid(spin_system,@coloc,parameters,'nmr')`.
- Lines 42-43: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'cos'},{'cos'}})`.
- Lines 45-47: Fourier transform; implemented by `spec=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 48-49: Plot; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=11.7`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'1H','13C'}`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={4.0 75.0}`.
- Lines 19: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=5.0`.
- Lines 20: computes `inter.coupling.scalar{2,2}` using `inter.coupling.scalar{2,2}=0`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `parameters.spins` using `parameters.spins={'1H','13C'}`.
- Lines 28: computes `parameters.delta2` using `parameters.delta2=30e-3`.
- Lines 29: computes `parameters.offset` using `parameters.offset=[2250 5000]`.
- Lines 30: computes `parameters.sweep` using `parameters.sweep=[5000 12000]`.
- Lines 31: computes `parameters.npoints` using `parameters.npoints=[256 256]`.
- Lines 32: computes `parameters.zerofill` using `parameters.zerofill=[512 512]`.
- Lines 33: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `fid` using `fid=liquid(spin_system,@coloc,parameters,'nmr')`.
- Lines 46-47: computes `spec` using `spec=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.

## Implementation structure

- A simple COLOC pulse sequence example for a two-spin
- 1H-13C system with a long-range J-coupling.
- Calculation time: seconds.
- Magnet field
- Spin system
- Interactions
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- Fourier transform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
