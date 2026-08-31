# examples/nmr_solids/rframe_nqi_fplanck.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/rframe_nqi_fplanck.m`
- Signature: `rframe_nqi_fplanck()`
- Total lines: 57

## Purpose

Powder magic angle spinning spectrum (rotor-synchronized detection) of a single quadrupolar 14N nucleus using 1D Fokker-Planck equation and a spherical grid. The calculation accounts for the second-order quadrupolar shift and lineshape by applying numerical second order corrections to the rotating frame transformation. Calculation time: hours

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Propagation is accelerated with a Krylov-subspace method, replacing direct matrix exponentiation by projection into a much smaller Arnoldi/Lanczos-type subspace.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.
- A Krylov-subspace or Arnoldi construction is used to avoid forming or exponentiating very large dense propagators directly.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: System specification; implemented by `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 17-18: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 21-22: Algorithmic options; implemented by `sys.disable={'trajlevel','krylov'}`.
- Lines 24-25: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: Experiment setup; implemented by `parameters.rate=50000`.
- Lines 44-45: Simulation; implemented by `fid=singlerot(spin_system,@acquire,parameters,'labframe')`.
- Lines 47-48: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 50-51: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 53-54: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 15: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(3.06e6,0.40,1,[0 0 0])`.
- Lines 18: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 19: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 22: computes `sys.disable` using `sys.disable={'trajlevel','krylov'}`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 29: computes `parameters.rate` using `parameters.rate=50000`.
- Lines 30: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 31: computes `parameters.max_rank` using `parameters.max_rank=85`.
- Lines 32: computes `parameters.grid` using `parameters.grid='rep_2ang_6400pts_sph'`.
- Lines 33: computes `parameters.sweep` using `parameters.sweep=50000`.
- Lines 34: computes `parameters.npoints` using `parameters.npoints=256`.
- Lines 35: computes `parameters.zerofill` using `parameters.zerofill=1024`.
- Lines 36: computes `parameters.offset` using `parameters.offset=18000`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'14N'}`.
- Lines 38: computes `parameters.rframes` using `parameters.rframes={{'14N',2}}`.
- Lines 39: computes `parameters.axis_units` using `parameters.axis_units='Hz'`.
- Lines 40: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','14N')`.

## Implementation structure

- Powder magic angle spinning spectrum (rotor-synchronized detection)
- of a single quadrupolar 14N nucleus using 1D Fokker-Planck equation
- and a spherical grid. The calculation accounts for the second-order
- quadrupolar shift and lineshape by applying numerical second order
- corrections to the rotating frame transformation.
- Calculation time: hours
- System specification
- Basis set
- Algorithmic options
- Spinach housekeeping
- Experiment setup
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
