# examples/nmr_liquids/pa_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/pa_strychnine.m`
- Signature: `pa_strychnine()`
- Total lines: 61

## Purpose

1H NMR spectrum of strychnine including an accurate model of line widths via Redfield superoperator. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Read spin system properties; implemented by `[sys,inter]=strychnine({'1H'})`.
- Lines 13-14: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 16-17: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 22-23: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 28-29: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 48-49: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 51-52: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 54-55: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 57-58: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H'})`.
- Lines 14: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 17: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 18: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 19: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 20: computes `bas.space_level` using `bas.space_level=1`.
- Lines 23: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 24: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 25: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 26: computes `inter.tau_c` using `inter.tau_c={200e-12}`.
- Lines 29: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 30: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 38: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 39: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 40: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 41: computes `parameters.offset` using `parameters.offset=2800`.

## Implementation structure

- 1H NMR spectrum of strychnine including an accurate model of
- line widths via Redfield superoperator.
- Calculation time: seconds
- Read spin system properties
- Magnet field
- Basis set
- Relaxation theory parameters
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
