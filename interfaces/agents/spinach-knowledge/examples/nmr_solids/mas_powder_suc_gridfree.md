# examples/nmr_solids/mas_powder_suc_gridfree.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_suc_gridfree.m`
- Signature: `mas_powder_suc_gridfree()`
- Total lines: 66

## Purpose

13C MAS spectrum of sucrose powder (assuming decoupling of 1H), computed using the grid-free Fokker-Planck MAS formalism. Che- mical shielding tensors, J-couplings and coordinates are esti- mated with DFT. A polyadic representation of the evolution ge- nerator is used, further particulars here: Calculation time: hours on a Tesla V100 GPU, much longer on CPU

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-18: Spin system properties (PCM DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'), {{'C','13C'}},182.1,[])`.
- Lines 19-20: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 22-23: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 28-29: Algorithmic options; implemented by `sys.tols.inter_cutoff=5.0`.
- Lines 33-34: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 37-38: Experiment setup; implemented by `parameters.axis=[1 1 1]`.
- Lines 53-54: Run the simulation; implemented by `fid=gridfree(spin_system,@acquire,parameters,'nmr')`.
- Lines 56-57: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 59-60: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 62-63: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 17-18: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'), {{'C','13C'}},182.1,[])`.
- Lines 20: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='IK-0'`.
- Lines 25: computes `bas.projections` using `bas.projections=+1`.
- Lines 26: computes `bas.level` using `bas.level=3`.
- Lines 29: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=5.0`.
- Lines 30: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 31: computes `sys.enable` using `sys.enable={'greedy','polyadic'}`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 38: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 39: computes `parameters.max_rank` using `parameters.max_rank=23`.
- Lines 40: computes `parameters.rate` using `parameters.rate=6000`.
- Lines 41: computes `parameters.sweep` using `parameters.sweep=5e4`.
- Lines 42: computes `parameters.npoints` using `parameters.npoints=256`.
- Lines 43: computes `parameters.zerofill` using `parameters.zerofill=1024`.
- Lines 44: computes `parameters.offset` using `parameters.offset=15000`.
- Lines 45: computes `parameters.spins` using `parameters.spins={'13C'}`.

## Implementation structure

- 13C MAS spectrum of sucrose powder (assuming decoupling of 1H),
- computed using the grid-free Fokker-Planck MAS formalism. Che-
- mical shielding tensors, J-couplings and coordinates are esti-
- mated with DFT. A polyadic representation of the evolution ge-
- nerator is used, further particulars here:
- Calculation time: hours on a Tesla V100 GPU,
- much longer on CPU
- Spin system properties (PCM DFT calculation)
- Magnet field
- Basis set
- Algorithmic options
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `gridfree()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
