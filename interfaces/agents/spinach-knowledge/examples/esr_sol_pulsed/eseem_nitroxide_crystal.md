# examples/esr_sol_pulsed/eseem_nitroxide_crystal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/eseem_nitroxide_crystal.m`
- Signature: `eseem_nitroxide_crystal()`
- Total lines: 68

## Purpose

Two-pulse X-band ESEEM spectrum of a nitroxide radical at a specific orientation relative to the lab frame. Magnetic parameters taken from a DFT calculation. Ideal pulses are assumed. Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Isotopes; implemented by `sys.isotopes={'E','14N'}`.
- Lines 15-16: Interactions; implemented by `inter.zeeman.matrix=cell(1,2)`.
- Lines 25-26: Magnet field; implemented by `sys.magnet=0.33`.
- Lines 28-29: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 47-48: Simulation; implemented by `fid=crystal(spin_system,@eseem,parameters,'esr')`.
- Lines 50-51: Plot the time domain signal; implemented by `kfigure(); subplot(2,1,1)`.
- Lines 55-56: Crude apodization; implemented by `fid=apodisation(spin_system,fid-mean(fid),{{'kaiser',6}})`.
- Lines 58-59: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 61-62: Plot the spectrum; implemented by `subplot(2,1,2)`.

### Key state/data transformations

- Lines 13: computes `sys.isotopes` using `sys.isotopes={'E','14N'}`.
- Lines 16: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,2)`.
- Lines 17: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=[2.01045 0.00000 0.00000`.
- Lines 20: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(2,2)`.
- Lines 21: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=[1.2356 0.0000 0.6322`.
- Lines 26: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 38: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 39: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 40: computes `parameters.screen` using `parameters.screen=state(spin_system,'L-','E')`.
- Lines 41: computes `parameters.pulse_op` using `parameters.pulse_op=operator(spin_system,'Ly','E')`.
- Lines 42: computes `parameters.npoints` using `parameters.npoints=1024`.
- Lines 43: computes `parameters.timestep` using `parameters.timestep=1.25e-8`.
- Lines 44: computes `parameters.orientation` using `parameters.orientation=[pi/5 pi/4 pi/3]`.
- Lines 45: computes `parameters.zerofill` using `parameters.zerofill=4096`.

## Implementation structure

- Two-pulse X-band ESEEM spectrum of a nitroxide radical at a specific
- orientation relative to the lab frame. Magnetic parameters taken from
- a DFT calculation. Ideal pulses are assumed.
- Calculation time: seconds
- Isotopes
- Interactions
- Magnet field
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Plot the time domain signal

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `crystal()`, `kfigure()`, `subplot()`, `kxlabel()`, `apodisation()`, `fftshift()`.
