# examples/nmr_solids/mas_powder_trp_gridfree.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_trp_gridfree.m`
- Signature: `mas_powder_trp_gridfree()`
- Total lines: 118

## Purpose

13C MAS spectrum of tryptophan powder (assuming decoupling of 1H), computed using the grid-free Fokker-Planck MAS formalism. Isotro- pic chemical shifts come from the experimental data. Coordinates and CSAs are estimated with DFT. A polyadic representation of the evolution generator is used, further particulars here: Calculation time: hours on a Tesla V100 GPU, much longer on CPU

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-22: % First molecule in the unit cell; implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/trp_xray.out'), {{'C','13C'},{'N','15N'}},[174.0 0],[])`.
- Lines 20-22: Spin system properties (DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/trp_xray.out'), {{'C','13C'},{'N','15N'}},[174.0 0],[])`.
- Lines 23-24: Magnet field; implemented by `sys.magnet=9.4`.
- Lines 26-27: First conformation; implemented by `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,2,124.2)`.
- Lines 39-40: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 46-47: Algorithmic options; implemented by `sys.tols.inter_cutoff=5.0`.
- Lines 51-52: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 55-56: Experiment setup; implemented by `parameters.rate=14000`.
- Lines 72-73: Run the simulation; implemented by `fid=gridfree(spin_system,@acquire,parameters,'nmr')`.
- Lines 75-79: % Second molecule in the unit cell; implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/trp_xray.out'), {{'C','13C'},{'N','15N'}},[174.0 0],[])`.
- Lines 88-89: Second conformation; implemented by `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,2,124.2)`.
- Lines 105-106: Run the simulation; implemented by `fid=fid+gridfree(spin_system,@acquire,parameters,'nmr')`.
- Lines 108-109: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 111-112: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 114-115: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 21-22: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/trp_xray.out'), {{'C','13C'},{'N','15N'}},[174.0 0],[])`.
- Lines 24: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 27: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,2,124.2)`.
- Lines 40: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 41: computes `bas.approximation` using `bas.approximation='IK-0'`.
- Lines 42: computes `bas.longitudinals` using `bas.longitudinals={'15N'}`.
- Lines 43: computes `bas.projections` using `bas.projections=+1`.
- Lines 44: computes `bas.level` using `bas.level=3`.
- Lines 47: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=5.0`.
- Lines 48: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 49: computes `sys.enable` using `sys.enable={'greedy','polyadic'}`.
- Lines 52: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 56: computes `parameters.rate` using `parameters.rate=14000`.
- Lines 57: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 58: computes `parameters.assumptions` using `parameters.assumptions='nmr'`.
- Lines 59: computes `parameters.max_rank` using `parameters.max_rank=11`.
- Lines 60: computes `parameters.sweep` using `parameters.sweep=1e5`.
- Lines 61: computes `parameters.npoints` using `parameters.npoints=2048`.

## Implementation structure

- 13C MAS spectrum of tryptophan powder (assuming decoupling of 1H),
- computed using the grid-free Fokker-Planck MAS formalism. Isotro-
- pic chemical shifts come from the experimental data. Coordinates
- and CSAs are estimated with DFT. A polyadic representation of the
- evolution generator is used, further particulars here:
- Calculation time: hours on a Tesla V100 GPU,
- much longer on CPU
- % First molecule in the unit cell
- Spin system properties (DFT calculation)
- Magnet field
- First conformation
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `basis()`, `state()`, `gridfree()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
