# examples/nmr_solids/mas_powder_dip_gridfree.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_dip_gridfree.m`
- Signature: `mas_powder_dip_gridfree()`
- Total lines: 54

## Purpose

Powder magic angle spinning spectrum of a pair of dipole-coupled proton spins using grid-free Fokker-Planck equation formalism. Calculation time: minutes

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: System specification; implemented by `sys.magnet=14.1`.
- Lines 16-17: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 21-22: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 25-26: Parameters; implemented by `parameters.axis=[1 1 1]`.
- Lines 41-42: Simulation; implemented by `fid=gridfree(spin_system,@acquire,parameters,'nmr')`.
- Lines 44-45: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 47-48: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 50-51: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 13: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={5.0,-2.0}`.
- Lines 14: computes `inter.coordinates` using `inter.coordinates={[0 0 0]; [0 3.9 0.1]}`.
- Lines 17: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 18: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 19: computes `bas.projections` using `bas.projections=+1`.
- Lines 22: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 26: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 27: computes `parameters.max_rank` using `parameters.max_rank=15`.
- Lines 28: computes `parameters.sweep` using `parameters.sweep=2e4`.
- Lines 29: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 30: computes `parameters.zerofill` using `parameters.zerofill=4096`.
- Lines 31: computes `parameters.offset` using `parameters.offset=0`.
- Lines 32: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 33: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 34: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 35: computes `parameters.invert_axis` using `parameters.invert_axis=1`.

## Implementation structure

- Powder magic angle spinning spectrum of a pair of dipole-coupled
- proton spins using grid-free Fokker-Planck equation formalism.
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
