# examples/nmr_solids/mas_powder_ala_gridfree.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_ala_gridfree.m`
- Signature: `mas_powder_ala_gridfree()`
- Total lines: 63

## Purpose

13C MAS spectrum of alanine powder (assuming decoupling of 1H), computed using the grid-free Fokker-Planck MAS formalism. All magnetic parameters are estimated from a DFT calculation. Calculation time: seconds on a Tesla V100 GPU, much longer on CPU

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-14: Spin system properties (PCM DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/alanine.log'), {{'C','13C'},{'N','15N'}},[182.1 264.5],[])`.
- Lines 15-16: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 18-19: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Algorithmic options; implemented by `sys.tols.inter_cutoff=5.0`.
- Lines 28-31: sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Experiment setup; implemented by `parameters.axis=[1 1 1]`.
- Lines 50-51: Simulation; implemented by `fid=gridfree(spin_system,@acquire,parameters,'nmr')`.
- Lines 53-54: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 56-57: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 59-60: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 13-14: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/alanine.log'), {{'C','13C'},{'N','15N'}},[182.1 264.5],[])`.
- Lines 16: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 19: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 20: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 21: computes `bas.longitudinals` using `bas.longitudinals={'15N'}`.
- Lines 22: computes `bas.projections` using `bas.projections=+1`.
- Lines 25: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=5.0`.
- Lines 26: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 27: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 36: computes `parameters.max_rank` using `parameters.max_rank=17`.
- Lines 37: computes `parameters.rate` using `parameters.rate=2000`.
- Lines 38: computes `parameters.sweep` using `parameters.sweep=5e4`.
- Lines 39: computes `parameters.npoints` using `parameters.npoints=256`.
- Lines 40: computes `parameters.zerofill` using `parameters.zerofill=1024`.
- Lines 41: computes `parameters.offset` using `parameters.offset=15000`.
- Lines 42: computes `parameters.spins` using `parameters.spins={'13C'}`.

## Implementation structure

- 13C MAS spectrum of alanine powder (assuming decoupling of 1H),
- computed using the grid-free Fokker-Planck MAS formalism. All
- magnetic parameters are estimated from a DFT calculation.
- Calculation time: seconds on a Tesla V100 GPU,
- much longer on CPU
- Spin system properties (PCM DFT calculation)
- Magnet field
- Basis set
- Algorithmic options
- sys.enable={'gpu'};
- Spinach housekeeping
- Experiment setup

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `gridfree()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
