# examples/esr_liq_pulsed/data_import/gaussian_import_example.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/data_import/gaussian_import_example.m`
- Signature: `gaussian_import_example()`
- Total lines: 54

## Purpose

Methyl radical simulation, Gaussian import. The uncommon signal intensity pattern comes from g-HFC cross-correlation.

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: System properties (vacuum DFT calculation); implemented by `options.no_xyz=1`.
- Lines 13-14: Magnet induction; implemented by `sys.magnet=0.339`.
- Lines 16-17: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 20-21: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 26-27: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 30-31: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 42-43: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'esr')`.
- Lines 45-46: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 48-49: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 51-52: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 9: computes `options.no_xyz` using `options.no_xyz=1`.
- Lines 10-11: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('gaussian_methyl_radical.out'), {{'E','E'},{'H','1H'}},[0 0],options)`.
- Lines 14: computes `sys.magnet` using `sys.magnet=0.339`.
- Lines 17: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 18: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 21: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 22: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 23: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 24: computes `inter.tau_c` using `inter.tau_c={5e-10}`.
- Lines 27: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 31: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 32: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','E')`.
- Lines 33: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 34: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 35: computes `parameters.sweep` using `parameters.sweep=5e8`.
- Lines 36: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 37: computes `parameters.zerofill` using `parameters.zerofill=1024`.
- Lines 38: computes `parameters.axis_units` using `parameters.axis_units='mT'`.

## Implementation structure

- Methyl radical simulation, Gaussian import. The uncommon signal
- intensity pattern comes from g-HFC cross-correlation.
- System properties (vacuum DFT calculation)
- Magnet induction
- Basis set
- Relaxation theory
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
