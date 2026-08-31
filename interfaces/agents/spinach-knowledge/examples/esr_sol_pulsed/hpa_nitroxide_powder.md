# examples/esr_sol_pulsed/hpa_nitroxide_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hpa_nitroxide_powder.m`
- Signature: `hpa_nitroxide_powder()`
- Total lines: 70

## Purpose

Powder averaged pulse-acquire W-band Fourier ESR spectrum of nitroxide radical. An ideal pulse is assumed. Calculation time: seconds

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
- Lines 26-27: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Disable trajectory-level SSR algorithms; implemented by `sys.disable={'trajlevel'}`.
- Lines 39-40: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 43-44: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 57-58: Simulation; implemented by `fid=powder(spin_system,@acquire,parameters,'esr')`.
- Lines 60-61: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'crisp'}})`.
- Lines 63-64: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 66-67: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11: computes `sys.isotopes` using `sys.isotopes={'E','14N'}`.
- Lines 14: computes `sys.magnet` using `sys.magnet=3.5`.
- Lines 17: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,2)`.
- Lines 18: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=[2.01045 0.00000 0.00000`.
- Lines 21: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(2,2)`.
- Lines 22: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=[1.2356 0.0000 0.6322`.
- Lines 27: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 28: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 29: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 30: computes `inter.damp_rate` using `inter.damp_rate=5e7`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 40: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 44: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 45: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','E')`.
- Lines 46: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 47: computes `parameters.decouple` using `parameters.decouple={}`.

## Implementation structure

- Powder averaged pulse-acquire W-band Fourier ESR spectrum of
- nitroxide radical. An ideal pulse is assumed.
- Calculation time: seconds
- Isotopes
- Magnet field
- Interactions
- Relaxation theory
- Basis set
- Disable trajectory-level SSR algorithms
- Spinach housekeeping
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
