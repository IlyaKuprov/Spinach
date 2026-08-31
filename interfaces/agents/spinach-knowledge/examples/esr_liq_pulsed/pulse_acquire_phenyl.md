# examples/esr_liq_pulsed/pulse_acquire_phenyl.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/pulse_acquire_phenyl.m`
- Signature: `pulse_acquire_phenyl()`
- Total lines: 61

## Purpose

W-band pulse-acquire FFT ESR spectrum of phenyl radical. Simple fixed line width is used as a relaxation model. Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Ignore coordinate information (HFCs provided); implemented by `options.no_xyz=1`.
- Lines 13-15: Read the spin system properties (vacuum DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/phenyl.log'), {{'E','E'},{'H','1H'}},[0 0],options)`.
- Lines 16-17: Magnet induction; implemented by `sys.magnet=3.5`.
- Lines 19-20: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 31-32: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 35-36: Set the sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 48-49: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'esr')`.
- Lines 51-52: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 54-55: Perform Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 57-58: Plot the spectrum; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 11: computes `options.no_xyz` using `options.no_xyz=1`.
- Lines 14-15: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/phenyl.log'), {{'E','E'},{'H','1H'}},[0 0],options)`.
- Lines 17: computes `sys.magnet` using `sys.magnet=3.5`.
- Lines 20: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 21: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 22: computes `bas.longitudinals` using `bas.longitudinals={'1H'}`.
- Lines 23: computes `bas.projections` using `bas.projections=+1`.
- Lines 26: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 27: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 28: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 29: computes `inter.damp_rate` using `inter.damp_rate=1e7`.
- Lines 32: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 36: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 37: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','E')`.
- Lines 38: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 39: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 40: computes `parameters.offset` using `parameters.offset=0`.
- Lines 41: computes `parameters.sweep` using `parameters.sweep=2e8`.

## Implementation structure

- W-band pulse-acquire FFT ESR spectrum of phenyl radical. Simple
- fixed line width is used as a relaxation model.
- Calculation time: seconds
- Ignore coordinate information (HFCs provided)
- Read the spin system properties (vacuum DFT calculation)
- Magnet induction
- Basis set
- Relaxation theory
- Spinach housekeeping
- Set the sequence parameters
- Simulation
- Apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
