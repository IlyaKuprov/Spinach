# examples/nmr_solids/static_powder_suc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/static_powder_suc.m`
- Signature: `static_powder_suc()`
- Total lines: 61

## Purpose

13C NMR spectrum of static sucrose powder, protons are assumed to be decoupled. Calculation time: hours

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-12: Spin system properties (PCM DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'), {{'C','13C'}},182.1,[])`.
- Lines 13-14: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 16-17: Disable trajectory-level SSR algorithms; implemented by `sys.disable={'trajlevel'}`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Algorithmic options; implemented by `sys.tols.inter_cutoff=5.0`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Experiment setup; implemented by `parameters.sweep=5e4`.
- Lines 48-49: Simulation; implemented by `fid=powder(spin_system,@acquire,parameters,'nmr')`.
- Lines 51-52: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 54-55: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 57-58: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11-12: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/sucrose.log'), {{'C','13C'}},182.1,[])`.
- Lines 14: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 17: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='IK-0'`.
- Lines 22: computes `bas.projections` using `bas.projections=+1`.
- Lines 23: computes `bas.level` using `bas.level=3`.
- Lines 26: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=5.0`.
- Lines 27: computes `sys.tols.prox_cutoff` using `sys.tols.prox_cutoff=4.0`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `parameters.sweep` using `parameters.sweep=5e4`.
- Lines 36: computes `parameters.npoints` using `parameters.npoints=128`.
- Lines 37: computes `parameters.zerofill` using `parameters.zerofill=512`.
- Lines 38: computes `parameters.offset` using `parameters.offset=15000`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'13C'}`.
- Lines 40: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 41: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 42: computes `parameters.invert_axis` using `parameters.invert_axis=1`.

## Implementation structure

- 13C NMR spectrum of static sucrose powder, protons are assumed
- to be decoupled.
- Calculation time: hours
- Spin system properties (PCM DFT calculation)
- Magnet field
- Disable trajectory-level SSR algorithms
- Basis set
- Algorithmic options
- Spinach housekeeping
- Experiment setup
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
