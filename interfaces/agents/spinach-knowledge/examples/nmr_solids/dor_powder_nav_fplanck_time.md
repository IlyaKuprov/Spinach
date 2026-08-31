# examples/nmr_solids/dor_powder_nav_fplanck_time.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/dor_powder_nav_fplanck_time.m`
- Signature: `dor_powder_nav_fplanck_time()`
- Total lines: 66

## Purpose

Double angle spinning spectrum of N-acetylvaline 14N nucleus using 1D Fokker-Planck equation and a spherical grid. The cal- culation includes the second-order quadrupolar shift and the third-order lineshape. Time-domain detection. Note: slower spinning rates and larger NQIs require larger ranks and spherical grids. At the moment the spinning frequencies are set artificially too high to reduce the simulation time in t

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: System specification; implemented by `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 21-22: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 27-28: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Algorithmic options; implemented by `sys.disable={'trajlevel'}`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Experiment setup; implemented by `parameters.rate_outer=1e6`.
- Lines 56-57: Simulation; implemented by `fid=doublerot(spin_system,@acquire,parameters,'labframe')`.
- Lines 59-60: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 62-63: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 18: computes `sys.magnet` using `sys.magnet=14.1; sys.isotopes={'14N'}`.
- Lines 19: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(3.21e6,0.27,1,[0 0 0])`.
- Lines 22: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 23: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 24: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 25: computes `inter.damp_rate` using `inter.damp_rate=2e3`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `parameters.rate_outer` using `parameters.rate_outer=1e6`.
- Lines 40: computes `parameters.rate_inner` using `parameters.rate_inner=5e6`.
- Lines 41: computes `parameters.rank_outer` using `parameters.rank_outer=7`.
- Lines 42: computes `parameters.rank_inner` using `parameters.rank_inner=4`.
- Lines 43: computes `parameters.axis_outer` using `parameters.axis_outer=[sqrt(2/3) 0 sqrt(1/3)]`.
- Lines 44: computes `parameters.axis_inner` using `parameters.axis_inner=[sqrt(20-2*sqrt(30)) 0 sqrt(15+2*sqrt(30))]`.
- Lines 45: computes `parameters.grid` using `parameters.grid='rep_2ang_100pts_oct'`.
- Lines 46: computes `parameters.sweep` using `parameters.sweep=1e5`.

## Implementation structure

- Double angle spinning spectrum of N-acetylvaline 14N nucleus
- using 1D Fokker-Planck equation and a spherical grid. The cal-
- culation includes the second-order quadrupolar shift and the
- third-order lineshape. Time-domain detection.
- Note: slower spinning rates and larger NQIs require larger
- ranks and spherical grids. At the moment the spinning
- frequencies are set artificially too high to reduce
- the simulation time in this example.
- Calculation time: seconds
- System specification
- Relaxation theory
- Basis set

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `doublerot()`, `fftshift()`, `kfigure()`, `plot_1d()`.
