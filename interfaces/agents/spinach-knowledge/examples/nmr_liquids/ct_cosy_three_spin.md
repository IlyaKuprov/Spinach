# examples/nmr_liquids/ct_cosy_three_spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/ct_cosy_three_spin.m`
- Signature: `ct_cosy_three_spin()`
- Total lines: 55

## Purpose

CT-COSY of three spin system. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Spin system and interactions; implemented by `sys.magnet=14.1`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 27-28: Sequence parameters; implemented by `parameters.offset=2700`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Simulation; implemented by `fid=liquid(spin_system,@ct_cosy,parameters,'nmr')`.
- Lines 42-43: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sqcos'},{'sqcos'}})`.
- Lines 45-47: Fourier transform; implemented by `spectrum=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 49-50: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H'}`.
- Lines 13: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.70 4.10 6.50}`.
- Lines 14: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=10`.
- Lines 15: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=8`.
- Lines 16: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=4`.
- Lines 17: computes `inter.coupling.scalar{3,3}` using `inter.coupling.scalar{3,3}=0`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 25: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 28: computes `parameters.offset` using `parameters.offset=2700`.
- Lines 29: computes `parameters.sweep` using `parameters.sweep=[3500 3500]`.
- Lines 30: computes `parameters.npoints` using `parameters.npoints=[256 256]`.
- Lines 31: computes `parameters.zerofill` using `parameters.zerofill=[512 512]`.
- Lines 32: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 33: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.

## Implementation structure

- CT-COSY of three spin system.
- Calculation time: seconds
- Spin system and interactions
- Basis set
- Algorithmic options
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
