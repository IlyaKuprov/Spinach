# examples/nmr_liquids/pa_sucrose.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/pa_sucrose.m`
- Signature: `pa_sucrose()`
- Total lines: 61

## Purpose

1H NMR spectrum of sucrose (magnetic parameters read in from a DFT calculation), including Redfield relaxation superoperator. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Read the spin system properties (vacuum DFT calculation); implemented by `options.min_j=1.0`.
- Lines 14-15: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 17-18: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 29-30: Proximity cut-off; implemented by `sys.tols.prox_cutoff=4.0`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 48-49: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 51-52: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 54-55: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 57-58: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11: computes `options.min_j` using `options.min_j=1.0`.
- Lines 12-13: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'), {{'H','1H'}},31.8,options)`.
- Lines 15: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 18: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 19: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 20: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 21: computes `bas.space_level` using `bas.space_level=2`.
- Lines 24: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 25: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 27: computes `inter.tau_c` using `inter.tau_c={1e-9}`.
- Lines 30: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 38: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 39: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 40: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 41: computes `parameters.offset` using `parameters.offset=1800`.

## Implementation structure

- 1H NMR spectrum of sucrose (magnetic parameters read in from a DFT
- calculation), including Redfield relaxation superoperator.
- Calculation time: seconds
- Read the spin system properties (vacuum DFT calculation)
- Magnet field
- Basis set
- Relaxation theory parameters
- Proximity cut-off
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
