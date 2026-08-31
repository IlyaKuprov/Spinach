# examples/nmr_liquids/hoesy_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/hoesy_strychnine.m`
- Signature: `hoesy_strychnine()`
- Total lines: 88

## Purpose

13C{1H} HOESY spectrum of strychnine with natural content of 13C isotope. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- The implementation explicitly addresses performance engineering through parallel or GPU execution, which matters because Spinach operators can become extremely large after basis expansion or powder/spatial lifting.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Strychnine spin system; implemented by `[sys,inter]=strychnine({'1H','13C'})`.
- Lines 14-15: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 17-18: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 31-32: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Sequence parameters; implemented by `parameters.tmix=0.5`.
- Lines 50-51: Generate isotopomers; implemented by `subsystems=dilute(spin_system,'13C')`.
- Lines 53-55: Preallocate the answer; implemented by `spectrum=zeros(parameters.zerofill(2), parameters.zerofill(1),'like',1i)`.
- Lines 57-58: Loop over isotopomers; implemented by `parfor n=1:numel(subsystems)`.
- Lines 60-61: Build the basis; implemented by `subsystem=basis(subsystems{n},bas)`.
- Lines 63-64: Simulation; implemented by `fid=liquid(subsystem,@hoesy,parameters,'nmr')`.
- Lines 66-67: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 70-71: F2 Fourier transform; implemented by `f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.
- Lines 74-75: Form States signal; implemented by `f1_states=f1_cos-1i*f1_sin`.
- Lines 77-78: F1 Fourier transform; implemented by `spectrum=spectrum+fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 82-83: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Control flow inferred from the code

- Line 58: `parfor` loop over `n=1:numel(subsystems)`.

### Key state/data transformations

- Lines 12: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H','13C'})`.
- Lines 15: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 18: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 19: computes `bas.approximation` using `bas.approximation='IK-1'`.
- Lines 20: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 21: computes `bas.space_level` using `bas.space_level=3`.
- Lines 22: computes `bas.level` using `bas.level=4`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 26: computes `inter.equilibrium` using `inter.equilibrium='IME'`.
- Lines 27: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 28: computes `inter.tau_c` using `inter.tau_c={50e-12}`.
- Lines 29: computes `inter.temperature` using `inter.temperature=298`.
- Lines 32: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 33: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=5.0`.
- Lines 34: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=2.0`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `parameters.tmix` using `parameters.tmix=0.5`.
- Lines 41: computes `parameters.sweep` using `parameters.sweep=[6000 18000]`.

## Implementation structure

- 13C{1H} HOESY spectrum of strychnine with natural content
- of 13C isotope.
- Calculation time: minutes
- Strychnine spin system
- Magnet field
- Basis set
- Relaxation theory parameters
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Generate isotopomers
- Preallocate the answer

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `dilute()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
