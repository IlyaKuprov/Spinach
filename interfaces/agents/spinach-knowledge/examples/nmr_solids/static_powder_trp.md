# examples/nmr_solids/static_powder_trp.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/static_powder_trp.m`
- Signature: `static_powder_trp()`
- Total lines: 73

## Purpose

13C NMR spectrum of tryptophan powder. Isotropic chemical shifts come from the experimental data. Coordinates and CSAs are estima- ted with DFT. Protons are assumed to be decoupled. Calculation time: hours

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-13: Spin system properties (DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/trp_xray.out'), {{'C','13C'},{'N','15N'}},[174.0 0],[])`.
- Lines 14-15: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 17-18: Experimental chemical shifts; implemented by `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,2,124.2)`.
- Lines 30-31: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 37-38: Algorithmic options; implemented by `sys.tols.inter_cutoff=5.0`.
- Lines 42-43: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 46-47: Experiment setup; implemented by `parameters.sweep=6e4`.
- Lines 60-61: Simulation; implemented by `fid=powder(spin_system,@acquire,parameters,'nmr')`.
- Lines 63-64: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 66-67: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 69-70: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 12-13: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/trp_xray.out'), {{'C','13C'},{'N','15N'}},[174.0 0],[])`.
- Lines 15: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 18: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=shift_iso(inter.zeeman.matrix,2,124.2)`.
- Lines 31: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 32: computes `bas.approximation` using `bas.approximation='IK-0'`.
- Lines 33: computes `bas.longitudinals` using `bas.longitudinals={'15N'}`.
- Lines 34: computes `bas.projections` using `bas.projections=+1`.
- Lines 35: computes `bas.level` using `bas.level=3`.
- Lines 38: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=5.0`.
- Lines 39: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 40: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 43: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 47: computes `parameters.sweep` using `parameters.sweep=6e4`.
- Lines 48: computes `parameters.npoints` using `parameters.npoints=128`.
- Lines 49: computes `parameters.zerofill` using `parameters.zerofill=512`.
- Lines 50: computes `parameters.offset` using `parameters.offset=18000`.
- Lines 51: computes `parameters.spins` using `parameters.spins={'13C'}`.
- Lines 52: computes `parameters.decouple` using `parameters.decouple={}`.

## Implementation structure

- 13C NMR spectrum of tryptophan powder. Isotropic chemical shifts
- come from the experimental data. Coordinates and CSAs are estima-
- ted with DFT. Protons are assumed to be decoupled.
- Calculation time: hours
- Spin system properties (DFT calculation)
- Magnet field
- Experimental chemical shifts
- Basis set
- Algorithmic options
- Spinach housekeeping
- Experiment setup
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `shift_iso()`, `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
