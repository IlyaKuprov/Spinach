# examples/esr_liq_pulsed/relaxation_parafluorotoluene.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/relaxation_parafluorotoluene.m`
- Signature: `relaxation_parafluorotoluene()`
- Total lines: 65

## Purpose

X-band pulse-acquire FFT ESR spectrum of parafluorotoluene radical, simulated using explicit time-domain propagation including Redfield relaxation superoperator. Calculation time: minutes

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Ignore coordinate information (HFCs provided); implemented by `options.no_xyz=1`.
- Lines 14-16: Read the spin system (vacuum DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/parafluorotoluene.log'), {{'E','E'},{'H','1H'},{'F','19F'}},[0 0 0],options)`.
- Lines 18-19: Ignore small HFC anisotropies; implemented by `sys.tols.inter_cutoff=1e5`.
- Lines 21-22: Magnet field; implemented by `sys.magnet=0.33`.
- Lines 24-25: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 29-30: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 35-36: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 39-40: Set the sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 52-53: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'esr')`.
- Lines 55-56: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 58-59: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 61-62: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 12: computes `options.no_xyz` using `options.no_xyz=1`.
- Lines 15-16: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/parafluorotoluene.log'), {{'E','E'},{'H','1H'},{'F','19F'}},[0 0 0],options)`.
- Lines 19: computes `sys.tols.inter_cutoff` using `sys.tols.inter_cutoff=1e5`.
- Lines 22: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 25: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `bas.zero_quantum` using `bas.zero_quantum={[1 2 3 5]}`.
- Lines 30: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 31: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 32: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 33: computes `inter.tau_c` using `inter.tau_c={1e-10}`.
- Lines 36: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 40: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 41: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','E')`.
- Lines 42: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 43: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 44: computes `parameters.offset` using `parameters.offset=0`.
- Lines 45: computes `parameters.sweep` using `parameters.sweep=3e8`.

## Implementation structure

- X-band pulse-acquire FFT ESR spectrum of parafluorotoluene
- radical, simulated using explicit time-domain propagation
- including Redfield relaxation superoperator.
- Calculation time: minutes
- Ignore coordinate information (HFCs provided)
- Read the spin system (vacuum DFT calculation)
- Ignore small HFC anisotropies
- Magnet field
- Basis set
- Relaxation theory
- Spinach housekeeping
- Set the sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
