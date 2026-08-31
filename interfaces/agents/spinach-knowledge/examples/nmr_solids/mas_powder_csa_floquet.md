# examples/nmr_solids/mas_powder_csa_floquet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_csa_floquet.m`
- Signature: `mas_powder_csa_floquet()`
- Total lines: 56

## Purpose

Powder magic angle spinning spectrum of a single anisotropically shielded proton spin using a Floquet theory based formalism. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file relies on Floquet theory, where periodic time dependence is lifted into an enlarged block representation that converts time-periodic dynamics into a time-independent eigenproblem.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: System specification; implemented by `sys.magnet=14.1`.
- Lines 17-18: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 22-23: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 26-27: Experiment setup; implemented by `parameters.rate=500`.
- Lines 43-44: Simulation; implemented by `fid=floquet(spin_system,@acquire,parameters,'nmr')`.
- Lines 46-47: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 49-50: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 52-53: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 14: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[-2 -2 4]-5,[-1 -3 4]+5}`.
- Lines 15: computes `inter.zeeman.euler` using `inter.zeeman.euler={[0 0 0],[0 0 0]}`.
- Lines 18: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 19: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 20: computes `bas.projections` using `bas.projections=+1`.
- Lines 23: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 27: computes `parameters.rate` using `parameters.rate=500`.
- Lines 28: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 29: computes `parameters.max_rank` using `parameters.max_rank=17`.
- Lines 30: computes `parameters.grid` using `parameters.grid='leb_2ang_rank_17'`.
- Lines 31: computes `parameters.sweep` using `parameters.sweep=2e4`.
- Lines 32: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 33: computes `parameters.zerofill` using `parameters.zerofill=4096`.
- Lines 34: computes `parameters.offset` using `parameters.offset=0`.
- Lines 35: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 36: computes `parameters.decouple` using `parameters.decouple={}`.

## Implementation structure

- Powder magic angle spinning spectrum of a single anisotropically shielded
- proton spin using a Floquet theory based formalism.
- Calculation time: seconds
- System specification
- Basis set
- Spinach housekeeping
- Experiment setup
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `floquet()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
