# examples/nmr_solids/mas_powder_csa_gridfree.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_csa_gridfree.m`
- Signature: `mas_powder_csa_gridfree()`
- Total lines: 55

## Purpose

Powder magic angle spinning spectrum of a pair of anisotropically shielded proton spins using grid- free Fokker-Planck equation formalism. Calculation time: minutes

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: System specification; implemented by `sys.magnet=14.1`.
- Lines 17-18: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 22-23: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 26-27: Parameters; implemented by `parameters.axis=[1 1 1]`.
- Lines 42-43: Simulation; implemented by `fid=gridfree(spin_system,@acquire,parameters,'nmr')`.
- Lines 45-46: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 48-49: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 51-52: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 13: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 14: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[-2 -2 4]-5,[-1 -3 4]+5}`.
- Lines 15: computes `inter.zeeman.euler` using `inter.zeeman.euler={[0 0 0],[0 0 0]}`.
- Lines 18: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 19: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 20: computes `bas.projections` using `bas.projections=+1`.
- Lines 23: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 27: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 28: computes `parameters.sweep` using `parameters.sweep=2e4`.
- Lines 29: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 30: computes `parameters.zerofill` using `parameters.zerofill=4096`.
- Lines 31: computes `parameters.offset` using `parameters.offset=0`.
- Lines 32: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 33: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 34: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 35: computes `parameters.invert_axis` using `parameters.invert_axis=1`.
- Lines 36: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.

## Implementation structure

- Powder magic angle spinning spectrum of a pair of
- anisotropically shielded proton spins using grid-
- free Fokker-Planck equation formalism.
- Calculation time: minutes
- System specification
- Basis set
- Spinach housekeeping
- Parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `gridfree()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
