# examples/nmr_solids/mqmas_nqi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_solids/mqmas_nqi.m`
- Signature: `mqmas_nqi()`
- Total lines: 59

## Purpose

Rotor-synchronous MQMAS spectrum of a 87Rb compound, transmitter set to the isotropic chemical shift. Calculation time: minutes

## Physical / mathematical content

- Solid-state NMR examples. The key physics is anisotropic spin interactions under static or magic-angle-spinning conditions: chemical-shift anisotropy, dipolar coupling, quadrupolar coupling, cross-polarisation, and orientation averaging using Floquet, Fokker-Planck, or direct powder quadrature formalisms.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: System specification: just the NQI; implemented by `sys.magnet=9.4; sys.isotopes={'87Rb'}`.
- Lines 14-15: Formalism and basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 18-19: Algorithmic options; implemented by `sys.disable={'trajlevel'}`.
- Lines 21-22: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 25-26: Experiment setup; implemented by `parameters.rate=62.5e3`.
- Lines 43-44: Simulation; implemented by `fid=singlerot(spin_system,@mqmas,parameters,'labframe')`.
- Lines 46-47: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'sqcos'},{'sqcos'}})`.
- Lines 49-51: Fourier transform; implemented by `spectrum=fftshift(fft2(fid,parameters.zerofill(2), parameters.zerofill(1)))`.
- Lines 53-54: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=9.4; sys.isotopes={'87Rb'}`.
- Lines 12: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(5e6,0.50,3/2,[0 0 0])`.
- Lines 15: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 16: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 19: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 22: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 26: computes `parameters.rate` using `parameters.rate=62.5e3`.
- Lines 27: computes `parameters.axis` using `parameters.axis=[1 1 1]`.
- Lines 28: computes `parameters.max_rank` using `parameters.max_rank=7`.
- Lines 29: computes `parameters.grid` using `parameters.grid='rep_2ang_1600pts_sph'`.
- Lines 30: computes `parameters.npoints` using `parameters.npoints=[128 128]`.
- Lines 31: computes `parameters.zerofill` using `parameters.zerofill=[256 256]`.
- Lines 32: computes `parameters.sweep` using `parameters.sweep=parameters.rate`.
- Lines 33: computes `parameters.offset` using `parameters.offset=0`.
- Lines 34: computes `parameters.spins` using `parameters.spins={'87Rb'}`.
- Lines 35: computes `parameters.rframes` using `parameters.rframes={{'87Rb',2}}`.
- Lines 36: computes `parameters.mq_order` using `parameters.mq_order=3`.
- Lines 37: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.

## Implementation structure

- Rotor-synchronous MQMAS spectrum of a 87Rb compound,
- transmitter set to the isotropic chemical shift.
- Calculation time: minutes
- System specification: just the NQI
- Formalism and basis set
- Algorithmic options
- Spinach housekeeping
- Experiment setup
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
