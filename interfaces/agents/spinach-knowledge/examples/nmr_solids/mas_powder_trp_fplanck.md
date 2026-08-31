# examples/nmr_solids/mas_powder_trp_fplanck.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_trp_fplanck.m`
- Signature: `mas_powder_trp_fplanck()`
- Total lines: 106

## Purpose

13C MAS spectrum of tryptophan powder (assuming decoupling of 1H), computed using the Fokker-Planck MAS formalism. Isotropic chemical shifts come from the experimental data. Coordinates are from X-ray data and CSAs are estimated with DFT. Calculation time: days, hours with a Tesla A100 GPU.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-18: % First molecule in the unit cell; implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/trp_xray.out'), {{'C','13C'},{'N','15N'}},[174.0 0],[])`.
- Lines 16-18: Spin system properties (DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/trp_xray.out'), {{'C','13C'},{'N','15N'}},[174.0 0],[])`.
- Lines 19-20: Magnet field; implemented by `sys.magnet=9.4`.
- Lines 22-23: Experimental chemical shifts, first conformation; implemented by `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,2,124.2)`.
- Lines 35-36: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 42-43: Algorithmic options; implemented by `sys.tols.inter_cutoff=5.0`.
- Lines 45-48: sys.enable={'gpu'};; implemented by `spin_system=create(sys,inter)`.
- Lines 47-48: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 51-52: Experiment setup; implemented by `parameters.rate=14000`.
- Lines 64-65: Simulation; implemented by `fid=singlerot(spin_system,@acquire,parameters,'nmr')`.
- Lines 67-71: % Second molecule in the unit cell; implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/trp_xray.out'), {{'C','13C'},{'N','15N'}},[174.0 0],[])`.
- Lines 80-81: Experimental chemical shifts, second conformation; implemented by `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,2,124.2)`.
- Lines 93-94: Add to the previous simulation; implemented by `fid=fid+singlerot(spin_system,@acquire,parameters,'nmr')`.
- Lines 96-97: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 99-100: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 102-103: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 17-18: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/trp_xray.out'), {{'C','13C'},{'N','15N'}},[174.0 0],[])`.
- Lines 20: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 23: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,2,124.2)`.
- Lines 36: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 37: computes `bas.approximation` using `bas.approximation='IK-0'`.
- Lines 38: computes `bas.longitudinals` using `bas.longitudinals={'15N'}`.
- Lines 39: computes `bas.projections` using `bas.projections=+1`.
- Lines 40: computes `bas.level` using `bas.level=3`.
- Lines 43: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=5.0`.
- Lines 44: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 48: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 52: computes `parameters.rate` using `parameters.rate=14000`.
- Lines 53: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 54: computes `parameters.max_rank` using `parameters.max_rank=11`.
- Lines 55: computes `parameters.grid` using `parameters.grid='leb_2ang_rank_11'`.
- Lines 56: computes `parameters.sweep` using `parameters.sweep=1e5`.
- Lines 57: computes `parameters.npoints` using `parameters.npoints=2048`.
- Lines 58: computes `parameters.zerofill` using `parameters.zerofill=8192`.

## Implementation structure

- 13C MAS spectrum of tryptophan powder (assuming decoupling of 1H),
- computed using the Fokker-Planck MAS formalism. Isotropic chemical
- shifts come from the experimental data. Coordinates are from X-ray
- data and CSAs are estimated with DFT.
- Calculation time: days, hours with a Tesla A100 GPU.
- % First molecule in the unit cell
- Spin system properties (DFT calculation)
- Magnet field
- Experimental chemical shifts, first conformation
- Basis set
- Algorithmic options
- sys.enable={'gpu'};

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
