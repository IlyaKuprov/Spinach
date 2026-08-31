# examples/nmr_solids/static_powder_nqi_b.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/static_powder_nqi_b.m`
- Signature: `static_powder_nqi_b()`
- Total lines: 69

## Purpose

Static powder 79Br NMR spectrum of potassium bromide. At least 3 quadrupolar tensors are necessary to reproduce the experimen- tal shape, likely due to a distribution of electrostatic envi- ronments in the powder. Calculation time: seconds.

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnet field; implemented by `sys.magnet=9.3659`.
- Lines 17-18: Spin system; implemented by `sys.isotopes={'79Br','79Br','79Br'}`.
- Lines 20-21: Chemical shift, ppm; implemented by `inter.zeeman.scalar={60.0933 60.0933 60.0933}`.
- Lines 23-24: Quadrupolar coupling; implemented by `inter.coupling.matrix{1,1}=1e3*diag([13.7569 1.6424 -(13.7569 + 1.6424)])`.
- Lines 28-29: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 33-34: Algorithmic options; implemented by `sys.disable={'trajlevel'}`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Experiment parameters; implemented by `parameters.grid='icos_2ang_163842pts'`.
- Lines 55-56: Simulation; implemented by `fid=powder(spin_system,@acquire,parameters,'nmr')`.
- Lines 58-59: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 61-62: Fourier transform; implemented by `spec=real(fftshift(fft(fid,parameters.zerofill)))`.
- Lines 64-65: Plotting; implemented by `kfigure(); plot_1d(spin_system,spec,parameters)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=9.3659`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'79Br','79Br','79Br'}`.
- Lines 21: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={60.0933 60.0933 60.0933}`.
- Lines 24: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=1e3*diag([13.7569 1.6424 -(13.7569 + 1.6424)])`.
- Lines 25: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=1e3*diag([4.0779 4.5179 -( 4.0779 + 4.5179)])`.
- Lines 26: computes `inter.coupling.matrix{3,3}` using `inter.coupling.matrix{3,3}=1e3*diag([1.5885 0.9449 -( 1.5885 + 0.9449)])`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='IK-0'`.
- Lines 31: computes `bas.level` using `bas.level=1; bas.projections=+1`.
- Lines 34: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `parameters.grid` using `parameters.grid='icos_2ang_163842pts'`.
- Lines 42: computes `parameters.spins` using `parameters.spins={'79Br'}`.
- Lines 43: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 44: computes `parameters.axis_units` using `parameters.axis_units='Hz'`.
- Lines 45: computes `parameters.invert_axis` using `parameters.invert_axis=1`.
- Lines 46: computes `parameters.offset` using `parameters.offset=6034.96`.
- Lines 47: computes `parameters.sweep` using `parameters.sweep=1e5`.

## Implementation structure

- Static powder 79Br NMR spectrum of potassium bromide. At least
- 3 quadrupolar tensors are necessary to reproduce the experimen-
- tal shape, likely due to a distribution of electrostatic envi-
- ronments in the powder.
- Calculation time: seconds.
- Magnet field
- Spin system
- Chemical shift, ppm
- Quadrupolar coupling
- Basis set
- Algorithmic options
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`, `ylim()`.
