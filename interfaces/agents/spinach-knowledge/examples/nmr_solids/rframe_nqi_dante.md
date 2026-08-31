# examples/nmr_solids/rframe_nqi_dante.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/rframe_nqi_dante.m`
- Signature: `rframe_nqi_dante()`
- Total lines: 65

## Purpose

DANTE MAS spectrum of a single quadrupolar 14N nucleus using 1D Fokker-Planck equation and a spherical grid. The calculation accounts for the second-order quadrupolar shift and lineshape. Set to reproduce Figure 3d from Calculation time: minutes

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: System specification; implemented by `sys.magnet=18.8; sys.isotopes={'14N'}`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Algorithmic options; implemented by `sys.disable={'trajlevel'}`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Experiment setup; implemented by `parameters.rate=62.5e3`.
- Lines 52-53: Simulation; implemented by `fid=singlerot(spin_system,@dante,parameters,'labframe')`.
- Lines 55-56: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 58-59: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 61-62: Plotting; implemented by `kfigure(); plot_1d(spin_system,abs(spectrum),parameters)`.

### Key state/data transformations

- Lines 17: computes `sys.magnet` using `sys.magnet=18.8; sys.isotopes={'14N'}`.
- Lines 18: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(1.18e6,0.50,1,[0 0 0])`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 26: computes `sys.enable` using `sys.enable={'prop_cache'}`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `parameters.rate` using `parameters.rate=62.5e3`.
- Lines 34: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 35: computes `parameters.max_rank` using `parameters.max_rank=35`.
- Lines 36: computes `parameters.grid` using `parameters.grid='rep_2ang_200pts_sph'`.
- Lines 37: computes `parameters.sweep` using `parameters.sweep=2000000`.
- Lines 38: computes `parameters.npoints` using `parameters.npoints=1024`.
- Lines 39: computes `parameters.zerofill` using `parameters.zerofill=4096`.
- Lines 40: computes `parameters.offset` using `parameters.offset=2200`.
- Lines 41: computes `parameters.spins` using `parameters.spins={'14N'}`.
- Lines 42: computes `parameters.rframes` using `parameters.rframes={{'14N',2}}`.
- Lines 43: computes `parameters.axis_units` using `parameters.axis_units='Hz'`.

## Implementation structure

- DANTE MAS spectrum of a single quadrupolar 14N nucleus using 1D
- Fokker-Planck equation and a spherical grid. The calculation
- accounts for the second-order quadrupolar shift and lineshape.
- Set to reproduce Figure 3d from
- Calculation time: minutes
- System specification
- Basis set
- Algorithmic options
- Spinach housekeeping
- Experiment setup
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
