# examples/esr_liq_pulsed/pulse_acquire_methyl.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/pulse_acquire_methyl.m`
- Signature: `pulse_acquire_methyl()`
- Total lines: 64

## Purpose

X-band pulse-acquire FFT ESR spectrum of methyl radical. Simple common line width is used as a relaxation model. Set to reprodu- ce Figure 4 from the paper by Zhitnikov and Dmitriev: Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Ignore coordinate information (HFCs provided); implemented by `options.no_xyz=1`.
- Lines 16-18: Read the spin system (vacuum DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/methyl.log'), {{'E','E'},{'H','1H'}},[0 0],options)`.
- Lines 19-20: Magnet induction; implemented by `sys.magnet=0.33`.
- Lines 22-23: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 28-29: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 34-35: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 38-39: Set the sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 51-52: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'esr')`.
- Lines 54-55: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 57-58: Perform Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 60-61: Plot the spectrum; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 14: computes `options.no_xyz` using `options.no_xyz=1`.
- Lines 17-18: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/methyl.log'), {{'E','E'},{'H','1H'}},[0 0],options)`.
- Lines 20: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `bas.projections` using `bas.projections=+1`.
- Lines 26: computes `bas.longitudinals` using `bas.longitudinals={'1H'}`.
- Lines 29: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 30: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 31: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 32: computes `inter.damp_rate` using `inter.damp_rate=2.5e7`.
- Lines 35: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 40: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','E')`.
- Lines 41: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 42: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 43: computes `parameters.offset` using `parameters.offset=0`.
- Lines 44: computes `parameters.sweep` using `parameters.sweep=5e8`.

## Implementation structure

- X-band pulse-acquire FFT ESR spectrum of methyl radical. Simple
- common line width is used as a relaxation model. Set to reprodu-
- ce Figure 4 from the paper by Zhitnikov and Dmitriev:
- Calculation time: seconds
- Ignore coordinate information (HFCs provided)
- Read the spin system (vacuum DFT calculation)
- Magnet induction
- Basis set
- Relaxation theory
- Spinach housekeeping
- Set the sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
