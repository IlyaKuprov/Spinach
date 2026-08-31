# examples/nmr_solids/cp_acquire_mas_gly.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_acquire_mas_gly.m`
- Signature: `cp_acquire_mas_gly()`
- Total lines: 80

## Purpose

1H-13C cross-polarisation followed by acquisition under magic angle spinning in alpha-glycine powder. Reduced Liouville spa- ce is used: up to, and including, three-spin correlations. Calculation time: minutes on Tesla A100, much longer on CPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-13: Spin system properties (PCM DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../../examples/standard_systems/glycine.log'), {{'H','1H'},{'C','13C'}},[31.8 182.1],[])`.
- Lines 15-16: 400 MHz spectrometer; implemented by `sys.magnet=9.4`.
- Lines 18-19: Isotropic alpha-glycine chemical shifts; implemented by `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1,176.4)`.
- Lines 27-28: Spin temperature; implemented by `inter.temperature=298`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 35-36: Algorithmic options; implemented by `sys.enable={'greedy'}`.
- Lines 39-40: Neglect interactions below 200 Hz; implemented by `sys.tols.inter_cutoff=2*pi*200`.
- Lines 42-43: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 46-47: Experiment setup; implemented by `parameters.spins={'1H','13C'}`.
- Lines 62-63: Detection state; implemented by `parameters.coil=state(spin_system,'L+','13C')`.
- Lines 65-66: Simulation; implemented by `fid=singlerot(spin_system,@cp_acquire_soft,parameters,'nmr')`.
- Lines 68-69: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 71-72: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 74-75: 1D plotting; implemented by `parameters.offset=parameters.offset(2)`.

### Key state/data transformations

- Lines 12-13: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../../examples/standard_systems/glycine.log'), {{'H','1H'},{'C','13C'}},[31.8 182.1],[])`.
- Lines 16: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 19: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,1,176.4)`.
- Lines 28: computes `inter.temperature` using `inter.temperature=298`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='IK-0'`.
- Lines 33: computes `bas.level` using `bas.level=3`.
- Lines 36: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 37: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 40: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=2*pi*200`.
- Lines 43: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 47: computes `parameters.spins` using `parameters.spins={'1H','13C'}`.
- Lines 48: computes `parameters.rate` using `parameters.rate=10000`.
- Lines 49: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 50: computes `parameters.max_rank` using `parameters.max_rank=5`.
- Lines 51: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 52: computes `parameters.offset` using `parameters.offset=[2e3 10e3]`.
- Lines 53: computes `parameters.hi_pwr` using `parameters.hi_pwr=83e3`.

## Implementation structure

- 1H-13C cross-polarisation followed by acquisition under magic
- angle spinning in alpha-glycine powder. Reduced Liouville spa-
- ce is used: up to, and including, three-spin correlations.
- Calculation time: minutes on Tesla A100, much longer on CPU.
- Spin system properties (PCM DFT calculation)
- 400 MHz spectrometer
- Isotropic alpha-glycine chemical shifts
- Spin temperature
- Basis set
- Algorithmic options
- Neglect interactions below 200 Hz
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
