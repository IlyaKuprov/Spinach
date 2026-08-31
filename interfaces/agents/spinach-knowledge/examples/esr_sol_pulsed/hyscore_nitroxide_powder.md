# examples/esr_sol_pulsed/hyscore_nitroxide_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hyscore_nitroxide_powder.m`
- Signature: `hyscore_nitroxide_powder()`
- Total lines: 65

## Purpose

Powder-averaged HYSCORE on a 14N nitroxide radical. Time-domain simulation in Liouville space. Set to reproduce Figure 2a from Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnet field; implemented by `sys.magnet=0.350`.
- Lines 15-16: System specification; implemented by `sys.isotopes={'14N','E'}`.
- Lines 23-24: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-28: Disable trajectory-level SSR algorithms; implemented by `sys.disable={'trajlevel','colorbar'}`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Set the sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 47-48: Simulation; implemented by `fid=powder(spin_system,@hyscore,parameters,'esr')`.
- Lines 50-51: Centre signal suppression; implemented by `fid=fid-mean(mean(fid))`.
- Lines 53-54: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'cos'},{'cos'}})`.
- Lines 56-57: Fourier transform; implemented by `spectrum=fftshift(fft2(fid,parameters.zerofill(2),parameters.zerofill(1)))`.
- Lines 59-60: Plotting; implemented by `kfigure(); scale_figure([1.5 2.0])`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=0.350`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'14N','E'}`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0 2.0000}`.
- Lines 18: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(2,2)`.
- Lines 19: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=5e6`.
- Lines 20: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(2,2)`.
- Lines 21: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=eeqq2nqi(2.4e6,0.5,1,[0 0 0])`.
- Lines 24: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 28: computes `sys.disable` using `sys.disable={'trajlevel','colorbar'}`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 36: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 37: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 38: computes `parameters.offset` using `parameters.offset=0`.
- Lines 39: computes `parameters.nsteps` using `parameters.nsteps=[128 128]`.
- Lines 40: computes `parameters.zerofill` using `parameters.zerofill=[256 256]`.
- Lines 41: computes `parameters.sweep` using `parameters.sweep=20e6`.

## Implementation structure

- Powder-averaged HYSCORE on a 14N nitroxide radical. Time-domain
- simulation in Liouville space. Set to reproduce Figure 2a from
- Calculation time: seconds
- Magnet field
- System specification
- Basis set
- Disable trajectory-level SSR algorithms
- Spinach housekeeping
- Set the sequence parameters
- Simulation
- Centre signal suppression
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `eeqq2nqi()`, `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `fft2()`, `kfigure()`, `scale_figure()`, `plot_2d()`.
