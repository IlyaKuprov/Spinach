# examples/nmr_solids/static_powder_ala.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/static_powder_ala.m`
- Signature: `static_powder_ala()`
- Total lines: 58

## Purpose

13C NMR spectrum of static alanine powder. Protons are assumed to be decoupled. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-12: Spin system properties (PCM DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/alanine.log'), {{'C','13C'},{'N','15N'}},[182.1 264.5],[])`.
- Lines 13-14: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 16-17: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 22-23: Algorithmic options; implemented by `sys.tols.inter_cutoff=5.0`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Experiment setup; implemented by `parameters.sweep=5e4`.
- Lines 45-46: Simulation; implemented by `fid=powder(spin_system,@acquire,parameters,'nmr')`.
- Lines 48-49: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 51-52: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 54-55: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11-12: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/alanine.log'), {{'C','13C'},{'N','15N'}},[182.1 264.5],[])`.
- Lines 14: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 17: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 18: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 19: computes `bas.longitudinals` using `bas.longitudinals={'15N'}`.
- Lines 20: computes `bas.projections` using `bas.projections=+1`.
- Lines 23: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=5.0`.
- Lines 24: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 25: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.sweep` using `parameters.sweep=5e4`.
- Lines 33: computes `parameters.npoints` using `parameters.npoints=128`.
- Lines 34: computes `parameters.zerofill` using `parameters.zerofill=512`.
- Lines 35: computes `parameters.offset` using `parameters.offset=15000`.
- Lines 36: computes `parameters.spins` using `parameters.spins={'13C'}`.
- Lines 37: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 38: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 39: computes `parameters.invert_axis` using `parameters.invert_axis=1`.

## Implementation structure

- 13C NMR spectrum of static alanine powder. Protons are assumed
- to be decoupled.
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

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
