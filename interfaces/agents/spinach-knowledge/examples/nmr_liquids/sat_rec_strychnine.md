# examples/nmr_liquids/sat_rec_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/sat_rec_strychnine.m`
- Signature: `sat_rec_strychnine()`
- Total lines: 62

## Purpose

1H saturation-recovery experiment on strychnine at 250 MHz. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Read spin system properties; implemented by `[sys,inter]=strychnine({'1H'})`.
- Lines 13-14: Magnetic induction; implemented by `sys.magnet=5.9`.
- Lines 16-17: Maximum distance to consider; implemented by `sys.tols.prox_cutoff=5.0`.
- Lines 19-20: Greedy parallelisation; implemented by `sys.enable={'greedy'}`.
- Lines 22-23: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 29-30: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 49-50: Simulation; implemented by `fids=liquid(spin_system,@sat_rec,parameters,'nmr')`.
- Lines 52-53: Apodisation; implemented by `fids=apodisation(spin_system,fids,{{'exp',6},{}})`.
- Lines 55-56: Fourier transform; implemented by `spectra=fftshift(fft(fids,[],1))`.
- Lines 58-59: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectra),parameters)`.

### Key state/data transformations

- Lines 11: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H'})`.
- Lines 14: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 17: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=5.0`.
- Lines 20: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 23: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 24: computes `inter.equilibrium` using `inter.equilibrium='dibari'`.
- Lines 25: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 26: computes `inter.tau_c` using `inter.tau_c={200e-12}`.
- Lines 27: computes `inter.temperature` using `inter.temperature=298`.
- Lines 30: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 31: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 32: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 33: computes `bas.space_level` using `bas.space_level=1`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 41: computes `parameters.offset` using `parameters.offset=1250`.
- Lines 42: computes `parameters.sweep` using `parameters.sweep=2500`.
- Lines 43: computes `parameters.npoints` using `parameters.npoints=4096`.

## Implementation structure

- 1H saturation-recovery experiment on strychnine at 250 MHz.
- Calculation time: minutes
- Read spin system properties
- Magnetic induction
- Maximum distance to consider
- Greedy parallelisation
- Relaxation theory parameters
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
