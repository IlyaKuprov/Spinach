# examples/nmr_solids/case_studies/magic_angle_calibration.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/case_studies/magic_angle_calibration.m`
- Signature: `magic_angle_calibration()`
- Total lines: 84

## Purpose

Magic angle is usually calibrated using KBr powder. When the angle is not correctly set, the spinning sideband pat- tern is blurred. This simulation demonstrates the effect. Calculation time: seconds

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnet field; implemented by `sys.magnet=9.4`.
- Lines 14-15: Isotopes; implemented by `sys.isotopes={'79Br'}`.
- Lines 17-18: Quardupolar coupling tensor; implemented by `inter.coupling.matrix{1,1}=eeqq2nqi(-92.4e3,-0.79,3/2,[0 0 0])`.
- Lines 20-21: Chemical shift; implemented by `inter.zeeman.scalar={60.0933}`.
- Lines 23-24: Formalism and basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Experiment setup; implemented by `parameters.spins={'79Br'}`.
- Lines 46-47: Convert magic angle errors to radians; implemented by `ma_errors=deg2rad([-1.00 -0.25 0.00 0.25 1.00])`.
- Lines 49-50: Get a figure going and compute the time axis; implemented by `kfigure(); scale_figure([2.5 2.5])`.
- Lines 53-54: Loop over angle errors; implemented by `for n=1:numel(ma_errors)`.
- Lines 56-57: Tilt the spinnig axis away from the magic angle; implemented by `parameters.axis=[sqrt(2/3) 0 sqrt(1/3)]*euler2dcm(0,ma_errors(n),0).'`.
- Lines 59-60: Run the MAS simulation; implemented by `fid=singlerot(spin_system,@acquire,parameters,'nmr')`.
- Lines 62-63: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 65-66: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 68-69: Plot the time domain signal; implemented by `subplot(numel(ma_errors),2,2*n-1)`.
- Lines 74-75: Plot the spectrum; implemented by `subplot(numel(ma_errors),2,2*n)`.

### Control flow inferred from the code

- Line 54: `for` loop over `n=1:numel(ma_errors)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'79Br'}`.
- Lines 18: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(-92.4e3,-0.79,3/2,[0 0 0])`.
- Lines 21: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={60.0933}`.
- Lines 24: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `bas.projections` using `bas.projections=+1`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `parameters.spins` using `parameters.spins={'79Br'}`.
- Lines 34: computes `parameters.rate` using `parameters.rate=4000`.
- Lines 35: computes `parameters.max_rank` using `parameters.max_rank=25`.
- Lines 36: computes `parameters.grid` using `parameters.grid='rep_2ang_1600pts_sph'`.
- Lines 37: computes `parameters.sweep` using `parameters.sweep=1e5`.
- Lines 38: computes `parameters.npoints` using `parameters.npoints=2048`.
- Lines 39: computes `parameters.zerofill` using `parameters.zerofill=4096`.
- Lines 40: computes `parameters.offset` using `parameters.offset=6000`.
- Lines 41: computes `parameters.axis_units` using `parameters.axis_units='Hz'`.
- Lines 42: computes `parameters.invert_axis` using `parameters.invert_axis=1`.

## Implementation structure

- Magic angle is usually calibrated using KBr powder. When
- the angle is not correctly set, the spinning sideband pat-
- tern is blurred. This simulation demonstrates the effect.
- Calculation time: seconds
- Magnet field
- Isotopes
- Quardupolar coupling tensor
- Chemical shift
- Formalism and basis set
- Spinach housekeeping
- Experiment setup
- Convert magic angle errors to radians

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `deg2rad()`, `kfigure()`, `scale_figure()`, `euler2dcm()`, `ma_errors()`, `singlerot()`, `apodisation()`, `fftshift()`, `subplot()`, `klegend()`, `num2str()`, `rad2deg()`.
