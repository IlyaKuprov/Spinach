# examples/nmr_liquids/roesy_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/roesy_strychnine.m`
- Signature: `roesy_strychnine()`
- Total lines: 69

## Purpose

ROESY spectrum of strychnine. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Spin system properties; implemented by `[sys,inter]=strychnine({'1H'})`.
- Lines 12-13: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 15-16: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 21-22: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 27-28: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Sequence parameters; implemented by `parameters.tmix=0.5`.
- Lines 46-47: Simulation; implemented by `fid=liquid(spin_system,@roesy,parameters,'nmr')`.
- Lines 49-50: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 53-54: F2 Fourier transform; implemented by `f1_cos=imag(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.
- Lines 57-58: States signal; implemented by `f1_states=f1_cos-1i*f1_sin`.
- Lines 60-61: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 63-64: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 10: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H'})`.
- Lines 13: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 16: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 17: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 18: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 19: computes `bas.space_level` using `bas.space_level=3`.
- Lines 22: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 23: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 24: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 25: computes `inter.tau_c` using `inter.tau_c={200e-12}`.
- Lines 28: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 29: computes `sys.disable` using `sys.disable={'krylov'}`.
- Lines 30: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.tmix` using `parameters.tmix=0.5`.
- Lines 38: computes `parameters.offset` using `parameters.offset=1200`.
- Lines 39: computes `parameters.sweep` using `parameters.sweep=[2500 2500]`.
- Lines 40: computes `parameters.npoints` using `parameters.npoints=[512 512]`.

## Implementation structure

- ROESY spectrum of strychnine.
- Calculation time: minutes
- Spin system properties
- Magnet field
- Basis set
- Relaxation theory parameters
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- F2 Fourier transform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
