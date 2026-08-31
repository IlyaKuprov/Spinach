# examples/esr_sol_pulsed/eseem_nitroxide_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/eseem_nitroxide_powder.m`
- Signature: `eseem_nitroxide_powder()`
- Total lines: 69

## Purpose

Powder-averaged two-pulse ESEEM on a 14N nitroxide radical. Time-domain simulation in Liouville space with powder averaging over a finite grid. Set to reproduce Figure 4a in http://dx.doi.org/10.1063/1.453532, ideal pulses are assumed. Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Magnet field; implemented by `sys.magnet=0.3249`.
- Lines 16-17: System specification; implemented by `sys.isotopes={'14N','E'}`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 29-30: Disable trajectory-level SSR algorithms; implemented by `sys.disable={'trajlevel'}`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Set the sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 48-49: Simulation; implemented by `fid=powder(spin_system,@eseem,parameters,'esr')`.
- Lines 51-52: Run apodization; implemented by `fid=apodisation(spin_system,mean(fid)-fid,{{'exp',5}})`.
- Lines 54-55: Run Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 57-58: Plot the time domain signal; implemented by `kfigure(); subplot(2,1,1)`.
- Lines 62-63: Plot the spectrum; implemented by `subplot(2,1,2)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=0.3249`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'14N','E'}`.
- Lines 18: computes `inter.coupling.eigs` using `inter.coupling.eigs=cell(2,2)`.
- Lines 19: computes `inter.coupling.euler` using `inter.coupling.euler=cell(2,2)`.
- Lines 20: computes `inter.coupling.eigs{1,1}` using `inter.coupling.eigs{1,1}=[-0.4 -1.6 +2.0]*1e5`.
- Lines 21: computes `inter.coupling.eigs{1,2}` using `inter.coupling.eigs{1,2}=[2.0 2.0 2.0]*1e6`.
- Lines 22: computes `inter.coupling.euler{1,1}` using `inter.coupling.euler{1,1}=[0 0 0]`.
- Lines 23: computes `inter.coupling.euler{1,2}` using `inter.coupling.euler{1,2}=[0 0 0]`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 30: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 38: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 39: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 40: computes `parameters.screen` using `parameters.screen=state(spin_system,'L-','E')`.
- Lines 41: computes `parameters.pulse_op` using `parameters.pulse_op=operator(spin_system,'Ly','E')`.
- Lines 42: computes `parameters.offset` using `parameters.offset=0`.

## Implementation structure

- Powder-averaged two-pulse ESEEM on a 14N nitroxide radical. Time-domain
- simulation in Liouville space with powder averaging over a finite grid.
- Set to reproduce Figure 4a in http://dx.doi.org/10.1063/1.453532, ideal
- pulses are assumed.
- Calculation time: seconds
- Magnet field
- System specification
- Basis set
- Disable trajectory-level SSR algorithms
- Spinach housekeeping
- Set the sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `subplot()`, `kxlabel()`.
