# examples/nmr_liquids/hoesy_ftyr_b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/hoesy_ftyr_b.m`
- Signature: `hoesy_ftyr_b()`
- Total lines: 74

## Purpose

(19F) -> (1H) HOESY spectrum of fluorotyrosine. This is not the right way to run this sequence in proteins because aro- matic 19F T2 is short, but 19F is being phase-encoded. Calculation time: minutes

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-13: Read 3-fluorotyrosine DFT calculation; implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/3_fluoro_tyr.log'), {{'H','1H'},{'F','19F'}},[31.82 192.97])`.
- Lines 15-16: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 18-19: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Relaxation theory parameters; implemented by `inter.relaxation={'redfield'}`.
- Lines 31-32: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Sequence parameters; implemented by `parameters.tmix=0.5`.
- Lines 51-52: Simulation; implemented by `fid=liquid(spin_system,@hoesy,parameters,'nmr')`.
- Lines 54-55: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 58-59: F2 Fourier transform; implemented by `f1_cos=real(fftshift(fft(fid.cos,parameters.zerofill(2),1),1))`.
- Lines 62-63: Form States signal; implemented by `f1_states=f1_cos-1i*f1_sin`.
- Lines 65-66: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 68-69: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 12-13: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/3_fluoro_tyr.log'), {{'H','1H'},{'F','19F'}},[31.82 192.97])`.
- Lines 16: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 19: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 20: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 21: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 22: computes `bas.space_level` using `bas.space_level=3`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 26: computes `inter.equilibrium` using `inter.equilibrium='IME'`.
- Lines 27: computes `inter.rlx_keep` using `inter.rlx_keep='kite'`.
- Lines 28: computes `inter.tau_c` using `inter.tau_c={10e-9}`.
- Lines 29: computes `inter.temperature` using `inter.temperature=298`.
- Lines 32: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 33: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=5.0`.
- Lines 34: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=2.0`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `parameters.tmix` using `parameters.tmix=0.5`.
- Lines 42: computes `parameters.sweep` using `parameters.sweep=[2500 4000]`.
- Lines 43: computes `parameters.offset` using `parameters.offset=[-70000 3000]`.

## Implementation structure

- (19F) -> (1H) HOESY spectrum of fluorotyrosine. This is not
- the right way to run this sequence in proteins because aro-
- matic 19F T2 is short, but 19F is being phase-encoded.
- Calculation time: minutes
- Read 3-fluorotyrosine DFT calculation
- Magnet field
- Basis set
- Relaxation theory parameters
- Algorithmic options
- Spinach housekeeping
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
