# examples/nmr_solids/mas_powder_dip_floquet.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_dip_floquet.m`
- Signature: `mas_powder_dip_floquet()`
- Total lines: 56

## Purpose

Spinning powder pulse-acquire experiment on a two-spin system with a dipolar coupling using Floquet theory. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file relies on Floquet theory, where periodic time dependence is lifted into an enlarged block representation that converts time-periodic dynamics into a time-independent eigenproblem.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: System specification; implemented by `sys.magnet=14.1`.
- Lines 17-18: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 22-23: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 26-27: Pulse-acquire setup; implemented by `parameters.rate=1000`.
- Lines 44-45: Simulation; implemented by `fid=floquet(spin_system,@acquire,parameters,'nmr')`.
- Lines 47-48: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 50-51: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 53-54: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 14: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={5.0,-2.0}`.
- Lines 15: computes `inter.coordinates` using `inter.coordinates={[0 0 0]; [0 3.9 0.1]}`.
- Lines 18: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 19: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 20: computes `bas.projections` using `bas.projections=+1`.
- Lines 23: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 27: computes `parameters.rate` using `parameters.rate=1000`.
- Lines 28: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 29: computes `parameters.max_rank` using `parameters.max_rank=17`.
- Lines 30: computes `parameters.grid` using `parameters.grid='leb_2ang_rank_17'`.
- Lines 31: computes `parameters.ref_frame` using `parameters.ref_frame='rotor'`.
- Lines 32: computes `parameters.sweep` using `parameters.sweep=2e4`.
- Lines 33: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 34: computes `parameters.zerofill` using `parameters.zerofill=4096`.
- Lines 35: computes `parameters.offset` using `parameters.offset=0`.
- Lines 36: computes `parameters.spins` using `parameters.spins={'1H'}`.

## Implementation structure

- Spinning powder pulse-acquire experiment on a two-spin system
- with a dipolar coupling using Floquet theory.
- Calculation time: seconds
- System specification
- Basis set
- Spinach housekeeping
- Pulse-acquire setup
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `floquet()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
