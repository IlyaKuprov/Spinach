# examples/esr_liq_pulsed/relaxation_nitroxide.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/relaxation_nitroxide.m`
- Signature: `relaxation_nitroxide()`
- Total lines: 60

## Purpose

W-band pulse-acquire FFT ESR spectrum of a nitroxide radical, using explicit time domain simulation with Redfield relaxati- on supeoperator. Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Ignore coordinate information (HFCs provided); implemented by `options.no_xyz=1`.
- Lines 14-16: Spin system properties (imported from a DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/nitroxide.log'), {{'E','E'},{'N','14N'}},[0 0],options)`.
- Lines 17-18: Magnet induction; implemented by `sys.magnet=3.5`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: RElaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 47-48: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'esr')`.
- Lines 50-51: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 53-54: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 56-57: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 12: computes `options.no_xyz` using `options.no_xyz=1`.
- Lines 15-16: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/nitroxide.log'), {{'E','E'},{'N','14N'}},[0 0],options)`.
- Lines 18: computes `sys.magnet` using `sys.magnet=3.5`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 26: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 27: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 28: computes `inter.tau_c` using `inter.tau_c={5e-11}`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 36: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','E')`.
- Lines 37: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 38: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 39: computes `parameters.offset` using `parameters.offset=-2e8`.
- Lines 40: computes `parameters.sweep` using `parameters.sweep=2e8`.
- Lines 41: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 42: computes `parameters.zerofill` using `parameters.zerofill=1024`.

## Implementation structure

- W-band pulse-acquire FFT ESR spectrum of a nitroxide radical,
- using explicit time domain simulation with Redfield relaxati-
- on supeoperator.
- Calculation time: seconds
- Ignore coordinate information (HFCs provided)
- Spin system properties (imported from a DFT calculation)
- Magnet induction
- Basis set
- RElaxation theory
- Spinach housekeeping
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
