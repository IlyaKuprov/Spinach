# examples/nmr_liquids/ct_cosy_2spins.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/ct_cosy_2spins.m`
- Signature: `ct_cosy_2spins()`
- Total lines: 53

## Purpose

CT COSY spectrum for 2 spins. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Spin system; implemented by `sys.isotopes={'1H','1H'}`.
- Lines 15-16: Interactions; implemented by `sys.magnet=5.9`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Sequence parameters; implemented by `parameters.offset=500`.
- Lines 33-34: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 37-38: Simulation; implemented by `fid=liquid(spin_system,@ct_cosy,parameters,'nmr')`.
- Lines 40-41: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sqcos'},{'sqcos'}})`.
- Lines 43-45: Fourier transform; implemented by `spectrum=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 47-48: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 16: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.00 5.00}`.
- Lines 18: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=7.0`.
- Lines 19: computes `inter.coupling.scalar{2,2}` using `inter.coupling.scalar{2,2}=0`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `parameters.offset` using `parameters.offset=500`.
- Lines 27: computes `parameters.sweep` using `parameters.sweep=[2000 2000]`.
- Lines 28: computes `parameters.npoints` using `parameters.npoints=[512 512]`.
- Lines 29: computes `parameters.zerofill` using `parameters.zerofill=[2048 2048]`.
- Lines 30: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 31: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 38: computes `fid` using `fid=liquid(spin_system,@ct_cosy,parameters,'nmr')`.
- Lines 44-45: computes `spectrum` using `spectrum=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.

## Implementation structure

- CT COSY spectrum for 2 spins.
- Calculation time: minutes
- Spin system
- Interactions
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
