# examples/nmr_liquids/cosy90_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/cosy90_strychnine.m`
- Signature: `cosy90_strychnine()`
- Total lines: 56

## Purpose

COSY spectrum of strychnine. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Read the spin system properties; implemented by `[sys,inter]=strychnine({'1H'})`.
- Lines 14-15: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 17-18: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Sequence parameters; implemented by `parameters.angle=pi/2`.
- Lines 40-41: Simulation; implemented by `fid=liquid(spin_system,@cosy,parameters,'nmr')`.
- Lines 43-44: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'cos'},{'cos'}})`.
- Lines 46-48: Fourier transform; implemented by `spectrum=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 50-51: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 12: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H'})`.
- Lines 15: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 18: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 19: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 24: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 25: computes `bas.space_level` using `bas.space_level=1`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.angle` using `parameters.angle=pi/2`.
- Lines 33: computes `parameters.offset` using `parameters.offset=1200`.
- Lines 34: computes `parameters.sweep` using `parameters.sweep=2200`.
- Lines 35: computes `parameters.npoints` using `parameters.npoints=[512 512]`.
- Lines 36: computes `parameters.zerofill` using `parameters.zerofill=[2048 2048]`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 38: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 41: computes `fid` using `fid=liquid(spin_system,@cosy,parameters,'nmr')`.
- Lines 47-48: computes `spectrum` using `spectrum=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.

## Implementation structure

- COSY spectrum of strychnine.
- Calculation time: minutes
- Read the spin system properties
- Magnet field
- Algorithmic options
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
