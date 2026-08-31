# examples/nmr_solids/mas_powder_dip_fplanck.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mas_powder_dip_fplanck.m`
- Signature: `mas_powder_dip_fplanck()`
- Total lines: 54

## Purpose

Spinning powder pulse-acquire experiment on a two-spin system with a dipolar coupling using Fokker-Planck formalism: Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: System specification; implemented by `sys.magnet=14.1`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: Pulse-acquire setup; implemented by `parameters.rate=1000`.
- Lines 41-42: Simulation; implemented by `fid=singlerot(spin_system,@acquire,parameters,'nmr')`.
- Lines 44-45: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 47-48: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 50-51: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={5.0,-2.0}`.
- Lines 16: computes `inter.coordinates` using `inter.coordinates={[0.0 0.0 0.0]`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 22: computes `bas.projections` using `bas.projections=+1`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 29: computes `parameters.rate` using `parameters.rate=1000`.
- Lines 30: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 31: computes `parameters.max_rank` using `parameters.max_rank=17`.
- Lines 32: computes `parameters.grid` using `parameters.grid='leb_2ang_rank_17'`.
- Lines 33: computes `parameters.sweep` using `parameters.sweep=2e4`.
- Lines 34: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 35: computes `parameters.zerofill` using `parameters.zerofill=4096`.
- Lines 36: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 37: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 38: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.

## Implementation structure

- Spinning powder pulse-acquire experiment on a two-spin system
- with a dipolar coupling using Fokker-Planck formalism:
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

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
