# examples/esr_sol_pulsed/eseem_methyl_crystal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/eseem_methyl_crystal.m`
- Signature: `eseem_methyl_crystal()`
- Total lines: 59

## Purpose

Two-pulse X-band ESEEM spectrum of a methyl radical at a specific orien- tation relative to the lab frame. Magnetic parameters taken from a DFT calculation. Ideal pulses are assumed. Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: System properties (vacuum DFT calculation); implemented by `options.no_xyz=1`.
- Lines 16-17: Magnet field; implemented by `sys.magnet=0.33`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 23-24: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 27-28: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 38-39: Simulation; implemented by `fid=crystal(spin_system,@eseem,parameters,'esr')`.
- Lines 41-42: Plot the time domain signal; implemented by `kfigure(); subplot(2,1,1)`.
- Lines 46-47: Crude apodization; implemented by `fid=apodisation(spin_system,fid-mean(fid),{{'kaiser',6}})`.
- Lines 49-50: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 52-53: Plot the spectrum; implemented by `subplot(2,1,2)`.

### Key state/data transformations

- Lines 13: computes `options.no_xyz` using `options.no_xyz=1`.
- Lines 14-15: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/methyl.log'), {{'E','E'},{'H','1H'}},[0 0],options)`.
- Lines 17: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 28: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 29: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 30: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 31: computes `parameters.screen` using `parameters.screen=state(spin_system,'L-','E')`.
- Lines 32: computes `parameters.pulse_op` using `parameters.pulse_op=operator(spin_system,'Ly','E')`.
- Lines 33: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 34: computes `parameters.timestep` using `parameters.timestep=1e-8`.
- Lines 35: computes `parameters.orientation` using `parameters.orientation=[pi/5 pi/4 pi/3]`.
- Lines 36: computes `parameters.zerofill` using `parameters.zerofill=4096`.
- Lines 39: computes `fid` using `fid=crystal(spin_system,@eseem,parameters,'esr')`.
- Lines 50: computes `spectrum` using `spectrum=fftshift(fft(fid,parameters.zerofill))`.

## Implementation structure

- Two-pulse X-band ESEEM spectrum of a methyl radical at a specific orien-
- tation relative to the lab frame. Magnetic parameters taken from a DFT
- calculation. Ideal pulses are assumed.
- Calculation time: seconds
- System properties (vacuum DFT calculation)
- Magnet field
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Plot the time domain signal
- Crude apodization

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `operator()`, `crystal()`, `kfigure()`, `subplot()`, `kxlabel()`, `apodisation()`, `fftshift()`.
