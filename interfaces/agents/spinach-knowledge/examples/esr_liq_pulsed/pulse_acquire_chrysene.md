# examples/esr_liq_pulsed/pulse_acquire_chrysene.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/pulse_acquire_chrysene.m`
- Signature: `pulse_acquire_chrysene()`
- Total lines: 68

## Purpose

W-band pulse-acquire FFT ESR spectrum of a chrysene cation radical in a non-viscous liquid. Simple common line width is used as a relaxation model. Symmetry treatment is performed using the full S2xS2xS2xS2xS2xS2 group direct product. Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Ignore coordinate information (HFCs provided); implemented by `options.no_xyz=1`.
- Lines 16-18: Read the spin system (vacuum DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/chrysene_cation.log'), {{'E','E'},{'H','1H'}},[0 0],options)`.
- Lines 19-20: Magnet induction; implemented by `sys.magnet=3.5`.
- Lines 22-23: Relaxation theory; implemented by `inter.relaxation={'damp'}`.
- Lines 28-29: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 34-35: Symmetry; implemented by `bas.sym_spins={[1 7],[2 8],[3 9],[4 10],[5 11],[6 12]}`.
- Lines 38-39: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 42-43: Set the sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 55-56: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'esr')`.
- Lines 58-59: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 61-62: Perform Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 64-65: Plot the spectrum; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 14: computes `options.no_xyz` using `options.no_xyz=1`.
- Lines 17-18: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/chrysene_cation.log'), {{'E','E'},{'H','1H'}},[0 0],options)`.
- Lines 20: computes `sys.magnet` using `sys.magnet=3.5`.
- Lines 23: computes `inter.relaxation` using `inter.relaxation={'damp'}`.
- Lines 24: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 25: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 26: computes `inter.damp_rate` using `inter.damp_rate=1e6`.
- Lines 29: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `bas.longitudinals` using `bas.longitudinals={'1H'}`.
- Lines 32: computes `bas.projections` using `bas.projections=+1`.
- Lines 35: computes `bas.sym_spins` using `bas.sym_spins={[1 7],[2 8],[3 9],[4 10],[5 11],[6 12]}`.
- Lines 36: computes `bas.sym_group` using `bas.sym_group={'S2','S2','S2','S2','S2','S2'}`.
- Lines 39: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 43: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 44: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','E')`.
- Lines 45: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 46: computes `parameters.decouple` using `parameters.decouple={}`.

## Implementation structure

- W-band pulse-acquire FFT ESR spectrum of a chrysene cation radical in
- a non-viscous liquid. Simple common line width is used as a relaxation
- model. Symmetry treatment is performed using the full S2xS2xS2xS2xS2xS2
- group direct product.
- Calculation time: seconds
- Ignore coordinate information (HFCs provided)
- Read the spin system (vacuum DFT calculation)
- Magnet induction
- Relaxation theory
- Basis set
- Symmetry
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
