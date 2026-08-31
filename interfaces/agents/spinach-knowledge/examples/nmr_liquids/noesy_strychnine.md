# examples/nmr_liquids/noesy_strychnine.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/noesy_strychnine.m`
- Signature: `noesy_strychnine()`
- Total lines: 71

## Purpose

NOESY spectrum of strychnine. Calculation time: minutes

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

- Lines 10-11: Spin system properties; implemented by `[sys,inter]=strychnine({'1H'})`.
- Lines 13-14: Magnet field; implemented by `sys.magnet=5.9`.
- Lines 16-17: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 22-23: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 29-30: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Sequence parameters; implemented by `parameters.tmix=0.5`.
- Lines 48-49: Simulation; implemented by `fid=liquid(spin_system,@noesy,parameters,'nmr')`.
- Lines 51-52: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 55-56: F2 Fourier transform; implemented by `f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.
- Lines 59-60: States signal; implemented by `f1_states=f1_cos-1i*f1_sin`.
- Lines 62-63: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 65-66: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 11: computes `[sys,inter]` using `[sys,inter]=strychnine({'1H'})`.
- Lines 14: computes `sys.magnet` using `sys.magnet=5.9`.
- Lines 17: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 18: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 19: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 20: computes `bas.space_level` using `bas.space_level=3`.
- Lines 23: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 24: computes `inter.equilibrium` using `inter.equilibrium='IME'`.
- Lines 25: computes `inter.temperature` using `inter.temperature=298`.
- Lines 26: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 27: computes `inter.tau_c` using `inter.tau_c={200e-12}`.
- Lines 30: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 31: computes `sys.disable` using `sys.disable={'krylov'}`.
- Lines 32: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `parameters.tmix` using `parameters.tmix=0.5`.
- Lines 40: computes `parameters.offset` using `parameters.offset=1200`.
- Lines 41: computes `parameters.sweep` using `parameters.sweep=[2500 2500]`.

## Implementation structure

- NOESY spectrum of strychnine.
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

- Called routines detected from the main body: `strychnine()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
