# examples/nmr_liquids/cosy90_derome.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/cosy90_derome.m`
- Signature: `cosy90_derome()`
- Total lines: 59

## Purpose

Figure 8.26 from Andrew Derome's "Modern NMR Techniques for Chemistry Research". Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnet field; implemented by `sys.magnet=16.1`.
- Lines 14-15: Spin system and interactions; implemented by `sys.isotopes={'1H','1H','1H'}`.
- Lines 22-23: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 26-27: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 30-31: Sequence parameters; implemented by `parameters.angle=pi/2`.
- Lines 39-40: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 43-44: Simulation; implemented by `fid=liquid(spin_system,@cosy,parameters,'nmr')`.
- Lines 46-47: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sqcos'},{'sqcos'}})`.
- Lines 49-51: Fourier transform; implemented by `spectrum=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 53-54: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=16.1`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H'}`.
- Lines 16: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={3.70 3.92 4.50}`.
- Lines 17: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=10`.
- Lines 18: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=12`.
- Lines 19: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=4`.
- Lines 20: computes `inter.coupling.scalar{3,3}` using `inter.coupling.scalar{3,3}=0`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 28: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 31: computes `parameters.angle` using `parameters.angle=pi/2`.
- Lines 32: computes `parameters.offset` using `parameters.offset=2800`.
- Lines 33: computes `parameters.sweep` using `parameters.sweep=700`.
- Lines 34: computes `parameters.npoints` using `parameters.npoints=[1024 1024]`.
- Lines 35: computes `parameters.zerofill` using `parameters.zerofill=[2048 2048]`.
- Lines 36: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 37: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.

## Implementation structure

- Figure 8.26 from Andrew Derome's "Modern NMR Techniques
- for Chemistry Research".
- Calculation time: seconds
- Magnet field
- Spin system and interactions
- Basis set
- Algorithmic options
- Sequence parameters
- Spinach housekeeping
- Simulation
- Apodisation
- Fourier transform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
