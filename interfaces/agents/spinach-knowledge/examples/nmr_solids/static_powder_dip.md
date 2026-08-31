# examples/nmr_solids/static_powder_dip.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/static_powder_dip.m`
- Signature: `static_powder_dip()`
- Total lines: 55

## Purpose

Static powder pulse-acquire experiment on a two-spin system with a dipolar coupling. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: System specification; implemented by `sys.magnet=14.1`.
- Lines 16-17: Algorithmic options; implemented by `sys.disable={'trajlevel'}`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 42-43: Simulation; implemented by `fid=powder(spin_system,@acquire,parameters,'nmr')`.
- Lines 45-46: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 48-49: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 51-52: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 13: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={5.0,-2.0}`.
- Lines 14: computes `inter.coordinates` using `inter.coordinates={[0 0 0]; [0 3.9 0.1]}`.
- Lines 17: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 22: computes `bas.projections` using `bas.projections=+1`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 29: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 30: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 31: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 32: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 33: computes `parameters.offset` using `parameters.offset=1000`.
- Lines 34: computes `parameters.grid` using `parameters.grid='rep_2ang_6400pts_sph'`.
- Lines 35: computes `parameters.sweep` using `parameters.sweep=12000`.
- Lines 36: computes `parameters.npoints` using `parameters.npoints=128`.
- Lines 37: computes `parameters.zerofill` using `parameters.zerofill=512`.

## Implementation structure

- Static powder pulse-acquire experiment on a two-spin system
- with a dipolar coupling.
- Calculation time: seconds
- System specification
- Algorithmic options
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
