# examples/nmr_solids/wise_mas_gly.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/wise_mas_gly.m`
- Signature: `wise_mas_gly()`
- Total lines: 82

## Purpose

WISE of alpha-glycine powder under MAS. Calculation time: hours, much faster on GPU

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-11: Spin system properties (PCM DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../../examples/standard_systems/glycine.log'), {{'H','1H'},{'C','13C'}},[31.8 182.1],[])`.
- Lines 13-14: 400 MHz spectrometer; implemented by `sys.magnet=9.4`.
- Lines 16-17: Isotropic alpha-glycine chemical shifts; implemented by `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1,176.4)`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Ignore interactions below 200 Hz; implemented by `sys.tols.inter_cutoff=2*pi*200`.
- Lines 33-34: Use GPU arithmetic; implemented by `sys.enable={'greedy'}`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Experiment setup; implemented by `parameters.spins={'1H','13C'}`.
- Lines 57-58: Detection state; implemented by `parameters.coil=state(spin_system,'L+','13C')`.
- Lines 60-61: Simulation; implemented by `fid=singlerot(spin_system,@wise,parameters,'nmr')`.
- Lines 63-64: Apodisation; implemented by `fid.cos=apodisation(spin_system,fid.cos,{{'sqcos'},{'sqcos'}})`.
- Lines 67-68: F2 Fourier transform; implemented by `f1_cos=fftshift(fft(fid.cos,parameters.zerofill(2),1),1)`.
- Lines 71-72: Form States signal; implemented by `f1_states=real(f1_cos)+1i*real(f1_sin)`.
- Lines 74-75: F1 Fourier transform; implemented by `spectrum=fftshift(fft(f1_states,parameters.zerofill(1),2),2)`.
- Lines 77-78: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 10-11: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../../examples/standard_systems/glycine.log'), {{'H','1H'},{'C','13C'}},[31.8 182.1],[])`.
- Lines 14: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 17: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1,176.4)`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='IK-0'`.
- Lines 28: computes `bas.level` using `bas.level=3`.
- Lines 31: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=2*pi*200`.
- Lines 34: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `parameters.spins` using `parameters.spins={'1H','13C'}`.
- Lines 42: computes `parameters.rate` using `parameters.rate=5000`.
- Lines 43: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 44: computes `parameters.max_rank` using `parameters.max_rank=9`.
- Lines 45: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 46: computes `parameters.offset` using `parameters.offset=[2e3 1e4]`.
- Lines 47: computes `parameters.hi_pwr` using `parameters.hi_pwr=83e3`.
- Lines 48: computes `parameters.cp_pwr` using `parameters.cp_pwr=[60e3 50e3]`.
- Lines 49: computes `parameters.cp_dur` using `parameters.cp_dur=1e-4`.

## Implementation structure

- WISE of alpha-glycine powder under MAS.
- Calculation time: hours, much faster on GPU
- Spin system properties (PCM DFT calculation)
- 400 MHz spectrometer
- Isotropic alpha-glycine chemical shifts
- Basis set
- Ignore interactions below 200 Hz
- Use GPU arithmetic
- Spinach housekeeping
- Experiment setup
- Detection state
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `stack_2d()`.
