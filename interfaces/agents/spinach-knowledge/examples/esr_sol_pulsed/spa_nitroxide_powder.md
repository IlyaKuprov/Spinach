# examples/esr_sol_pulsed/spa_nitroxide_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/spa_nitroxide_powder.m`
- Signature: `spa_nitroxide_powder()`
- Total lines: 73

## Purpose

A soft pulse simulation for a nitroxide radical powder. The soft pulse is simulated using the Fokker-Planck formalism; it is fol- lowed by time domain acquisition and Fourier transform. Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Isotopes; implemented by `sys.isotopes={'E','14N'}`.
- Lines 14-15: Magnet field; implemented by `sys.magnet=3.5`.
- Lines 17-18: Interactions; implemented by `inter.zeeman.matrix=cell(1,2)`.
- Lines 27-28: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Disable trajectory-level SSR algorithms; implemented by `sys.disable={'trajlevel'}`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 52-53: Soft pulse parameters; implemented by `parameters.pulse_rnk=2`.
- Lines 60-61: Simulation; implemented by `fid=powder(spin_system,@sp_acquire,parameters,'esr')`.
- Lines 63-64: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'crisp'}})`.
- Lines 66-67: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 69-70: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'E','14N'}`.
- Lines 15: computes `sys.magnet` using `sys.magnet=3.5`.
- Lines 18: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,2)`.
- Lines 19: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=[2.01045 0.00000 0.00000`.
- Lines 22: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(2,2)`.
- Lines 23: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=[1.2356 0.0000 0.6322`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 40: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 41: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 42: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 43: computes `parameters.offset` using `parameters.offset=-2e8`.
- Lines 44: computes `parameters.sweep` using `parameters.sweep=8e8`.
- Lines 45: computes `parameters.npoints` using `parameters.npoints=64`.
- Lines 46: computes `parameters.zerofill` using `parameters.zerofill=512`.

## Implementation structure

- A soft pulse simulation for a nitroxide radical powder. The soft
- pulse is simulated using the Fokker-Planck formalism; it is fol-
- lowed by time domain acquisition and Fourier transform.
- Calculation time: seconds
- Isotopes
- Magnet field
- Interactions
- Basis set
- Disable trajectory-level SSR algorithms
- Spinach housekeeping
- Sequence parameters
- Soft pulse parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
