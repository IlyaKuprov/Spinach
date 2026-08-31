# examples/nmr_solids/mas_powder_gly_floquet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_gly_floquet.m`
- Signature: `mas_powder_gly_floquet()`
- Total lines: 61

## Purpose

13C MAS spectrum of glycine powder (assuming decoupling of 1H), computed using the Floquet MAS formalism. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file relies on Floquet theory, where periodic time dependence is lifted into an enlarged block representation that converts time-periodic dynamics into a time-independent eigenproblem.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-13: Spin system properties (PCM DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/glycine.log'), {{'C','13C'},{'N','15N'}},[182.1 264.5],[])`.
- Lines 14-15: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 17-18: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Algorithmic options; implemented by `sys.tols.inter_cutoff=5.0`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Experiment setup; implemented by `parameters.rate=2000`.
- Lines 48-49: Simulation; implemented by `fid=floquet(spin_system,@acquire,parameters,'nmr')`.
- Lines 51-52: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 54-55: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 57-58: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 12-13: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/glycine.log'), {{'C','13C'},{'N','15N'}},[182.1 264.5],[])`.
- Lines 15: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 18: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 19: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 20: computes `bas.projections` using `bas.projections=+1`.
- Lines 21: computes `bas.longitudinals` using `bas.longitudinals={'15N'}`.
- Lines 24: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=5.0`.
- Lines 25: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.rate` using `parameters.rate=2000`.
- Lines 33: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 34: computes `parameters.max_rank` using `parameters.max_rank=23`.
- Lines 35: computes `parameters.sweep` using `parameters.sweep=5e4`.
- Lines 36: computes `parameters.npoints` using `parameters.npoints=256`.
- Lines 37: computes `parameters.zerofill` using `parameters.zerofill=1024`.
- Lines 38: computes `parameters.offset` using `parameters.offset=17000`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'13C'}`.
- Lines 40: computes `parameters.decouple` using `parameters.decouple={}`.

## Implementation structure

- 13C MAS spectrum of glycine powder (assuming decoupling of 1H),
- computed using the Floquet MAS formalism.
- Calculation time: seconds
- Spin system properties (PCM DFT calculation)
- Magnet field
- Basis set
- Algorithmic options
- Spinach housekeeping
- Experiment setup
- Simulation
- Apodisation
- Fourier transform

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `floquet()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
