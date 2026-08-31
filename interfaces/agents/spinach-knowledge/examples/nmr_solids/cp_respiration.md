# examples/nmr_solids/cp_respiration.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/cp_respiration.m`
- Signature: `cp_respiration()`
- Total lines: 57

## Purpose

1H-13C RESPIRATION-CP experiment in the doubly rotating frame. Magic angle spinning simulation using Fokker-Planck formalism. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnet field; implemented by `sys.magnet=11.7`.
- Lines 15-16: System specification; implemented by `sys.isotopes={'1H','13C'}`.
- Lines 20-21: Formalism and basis; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 28-29: Experiment parameters; implemented by `parameters.rate=20000`.
- Lines 43-44: Simulation; implemented by `fid=singlerot(spin_system,@respiration,parameters,'nmr')`.
- Lines 46-47: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',5}})`.
- Lines 49-50: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 52-53: Plotting; implemented by `kfigure(); parameters.spins={'13C'}`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=11.7`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'1H','13C'}`.
- Lines 17: computes `inter.coordinates` using `inter.coordinates={[0.00 0.00 0.00]`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 29: computes `parameters.rate` using `parameters.rate=20000`.
- Lines 30: computes `parameters.axis` using `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 31: computes `parameters.max_rank` using `parameters.max_rank=8`.
- Lines 32: computes `parameters.spins` using `parameters.spins={'1H','13C'}`.
- Lines 33: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lx','1H')`.
- Lines 34: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','13C')`.
- Lines 35: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_sph'`.
- Lines 36: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 37: computes `parameters.zerofill` using `parameters.zerofill=16384`.
- Lines 38: computes `parameters.axis_units` using `parameters.axis_units='kHz'`.
- Lines 39: computes `parameters.sweep` using `parameters.sweep=40000`.
- Lines 40: computes `parameters.nloops` using `parameters.nloops=16`.

## Implementation structure

- 1H-13C RESPIRATION-CP experiment in the doubly rotating frame.
- Magic angle spinning simulation using Fokker-Planck formalism.
- Calculation time: seconds
- Magnet field
- System specification
- Formalism and basis
- Spinach housekeeping
- Experiment parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
