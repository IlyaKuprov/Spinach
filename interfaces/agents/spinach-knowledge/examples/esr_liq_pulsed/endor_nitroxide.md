# examples/esr_liq_pulsed/endor_nitroxide.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/endor_nitroxide.m`
- Signature: `endor_nitroxide()`
- Total lines: 53

## Purpose

Mims ENDOR on a 15N-labelled nitroxide radical in liquid state. Magnetic parameters taken from a DFT calculation. Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Ignore coordinate information (HFCs provided); implemented by `options.no_xyz=1`.
- Lines 14-16: Read the spin system properties (vacuum DFT calculation); implemented by `[sys,inter]=g2spinach(gparse('../standard_systems/nitroxide.log'), {{'E','E'},{'N','15N'}},[0 0],options)`.
- Lines 17-18: Magnet field; implemented by `sys.magnet=0.33`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Disable path tracing (small system); implemented by `sys.disable={'pt'}`.
- Lines 27-28: Sequence parameters; implemented by `parameters.offset=0`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Simulation; implemented by `fid=liquid(spin_system,@endor_mims,parameters,'esr')`.
- Lines 43-44: Crude apodisation; implemented by `fid=apodisation(spin_system,fid-mean(fid),{{'kaiser',6}})`.
- Lines 46-47: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 49-50: Plotting; implemented by `kfigure(); plot_1d(spin_system,abs(spectrum),parameters)`.

### Key state/data transformations

- Lines 12: computes `options.no_xyz` using `options.no_xyz=1`.
- Lines 15-16: computes `[sys,inter]` using `[sys,inter]=g2spinach(gparse('../standard_systems/nitroxide.log'), {{'E','E'},{'N','15N'}},[0 0],options)`.
- Lines 18: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `sys.disable` using `sys.disable={'pt'}`.
- Lines 28: computes `parameters.offset` using `parameters.offset=0`.
- Lines 29: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 30: computes `parameters.sweep` using `parameters.sweep=1e8`.
- Lines 31: computes `parameters.tau` using `parameters.tau=100e-9`.
- Lines 32: computes `parameters.zerofill` using `parameters.zerofill=4096`.
- Lines 33: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 34: computes `parameters.axis_units` using `parameters.axis_units='MHz'`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `fid` using `fid=liquid(spin_system,@endor_mims,parameters,'esr')`.
- Lines 47: computes `spectrum` using `spectrum=fftshift(fft(fid,parameters.zerofill))`.

## Implementation structure

- Mims ENDOR on a 15N-labelled nitroxide radical in liquid
- state. Magnetic parameters taken from a DFT calculation.
- Calculation time: seconds
- Ignore coordinate information (HFCs provided)
- Read the spin system properties (vacuum DFT calculation)
- Magnet field
- Basis set
- Disable path tracing (small system)
- Sequence parameters
- Spinach housekeeping
- Simulation
- Crude apodisation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `g2spinach()`, `gparse()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
