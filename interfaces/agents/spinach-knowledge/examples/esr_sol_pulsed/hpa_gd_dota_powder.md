# examples/esr_sol_pulsed/hpa_gd_dota_powder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hpa_gd_dota_powder.m`
- Signature: `hpa_gd_dota_powder()`
- Total lines: 61

## Purpose

Powder averaged W-band pulsed ESR spectrum of Gd(III) DOTA complex. Ideal pulse with a large numerical powder grid is used, along with the numerical second-order rotating frame transformation. Calculation time: minutes

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Spin system properties; implemented by `sys.isotopes={'E8'}`.
- Lines 19-20: Magnet field; implemented by `sys.magnet=9.40`.
- Lines 22-23: Basis set; implemented by `bas.formalism='zeeman-liouv'`.
- Lines 26-27: Disable trajectory-level SSR algorithms; implemented by `sys.disable={'trajlevel'}`.
- Lines 29-30: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 33-34: Sequence parameters; implemented by `parameters.spins={'E8'}`.
- Lines 48-49: Simulation; implemented by `fid=powder(spin_system,@acquire,parameters,'labframe')`.
- Lines 51-52: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 54-55: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 57-58: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 14: computes `sys.isotopes` using `sys.isotopes={'E8'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.9918}`.
- Lines 16: computes `inter.coupling.eigs{1,1}` using `inter.coupling.eigs{1,1}=[0.57e9 0.57e9 -2*0.57e9]/3`.
- Lines 17: computes `inter.coupling.euler{1,1}` using `inter.coupling.euler{1,1}=[0 0 0]`.
- Lines 20: computes `sys.magnet` using `sys.magnet=9.40`.
- Lines 23: computes `bas.formalism` using `bas.formalism='zeeman-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 30: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 34: computes `parameters.spins` using `parameters.spins={'E8'}`.
- Lines 35: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','E8')`.
- Lines 36: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E8')`.
- Lines 37: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 38: computes `parameters.offset` using `parameters.offset=1.5e9`.
- Lines 39: computes `parameters.sweep` using `parameters.sweep=6e9`.
- Lines 40: computes `parameters.npoints` using `parameters.npoints=4096`.
- Lines 41: computes `parameters.zerofill` using `parameters.zerofill=16384`.
- Lines 42: computes `parameters.axis_units` using `parameters.axis_units='GHz-labframe'`.

## Implementation structure

- Powder averaged W-band pulsed ESR spectrum of Gd(III) DOTA
- complex. Ideal pulse with a large numerical powder grid is
- used, along with the numerical second-order rotating frame
- transformation.
- Calculation time: minutes
- Spin system properties
- Magnet field
- Basis set
- Disable trajectory-level SSR algorithms
- Spinach housekeeping
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `powder()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
