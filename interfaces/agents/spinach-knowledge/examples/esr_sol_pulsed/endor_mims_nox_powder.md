# examples/esr_sol_pulsed/endor_mims_nox_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/endor_mims_nox_powder.m`
- Signature: `endor_mims_nox_powder()`
- Total lines: 61

## Purpose

Mims ENDOR simulation for a nitroxide radical powder. Ideal hard pulses are assumed. Calculation time: seconds.

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Isotopes; implemented by `sys.isotopes={'E','14N'}`.
- Lines 13-14: Magnet field; implemented by `sys.magnet=3.5`.
- Lines 16-17: Interactions; implemented by `inter.zeeman.matrix=cell(1,2)`.
- Lines 26-27: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 30-31: Disable trajectory-level SSR algorithms; implemented by `sys.disable={'trajlevel'}`.
- Lines 33-34: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 37-38: Sequence parameters; implemented by `parameters.offset=0`.
- Lines 47-48: Simulation; implemented by `fid=powder(spin_system,@endor_mims,parameters,'esr')`.
- Lines 50-51: Crude apodisation; implemented by `fid=apodisation(spin_system,fid-mean(fid),{{'exp',6}})`.
- Lines 53-54: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 56-57: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11: computes `sys.isotopes` using `sys.isotopes={'E','14N'}`.
- Lines 14: computes `sys.magnet` using `sys.magnet=3.5`.
- Lines 17: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,2)`.
- Lines 18: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=[2.01045 0.00000 0.00000`.
- Lines 21: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(2,2)`.
- Lines 22: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=[1.2356 0.0000 0.6322`.
- Lines 27: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 38: computes `parameters.offset` using `parameters.offset=0`.
- Lines 39: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 40: computes `parameters.sweep` using `parameters.sweep=10e8`.
- Lines 41: computes `parameters.tau` using `parameters.tau=100e-9`.
- Lines 42: computes `parameters.zerofill` using `parameters.zerofill=4096`.
- Lines 43: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 44: computes `parameters.axis_units` using `parameters.axis_units='MHz'`.
- Lines 45: computes `parameters.grid` using `parameters.grid='rep_2ang_12800pts_sph'`.

## Implementation structure

- Mims ENDOR simulation for a nitroxide radical powder. Ideal
- hard pulses are assumed.
- Calculation time: seconds.
- Isotopes
- Magnet field
- Interactions
- Basis set
- Disable trajectory-level SSR algorithms
- Spinach housekeeping
- Sequence parameters
- Simulation
- Crude apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`, `kxlabel()`.
