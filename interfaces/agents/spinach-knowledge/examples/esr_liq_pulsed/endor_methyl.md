# examples/esr_liq_pulsed/endor_methyl.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/endor_methyl.m`
- Signature: `endor_methyl()`
- Total lines: 53

## Purpose

Mims ENDOR spectrum of a methyl radical in liquid state. Magnetic parameters taken from a DFT calculation. Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: Magnet field; implemented by `sys.magnet=0.346`.
- Lines 13-14: Spin system and interactions; implemented by `sys.isotopes={'1H','1H','1H','E'}`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-28: Sequence parameters; implemented by `parameters.offset=0`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Simulation; implemented by `fid=liquid(spin_system,@endor_mims,parameters,'esr')`.
- Lines 43-44: Crude apodisation; implemented by `fid=apodisation(spin_system,fid-mean(fid),{{'kaiser',6}})`.
- Lines 46-47: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 49-50: Plotting; implemented by `kfigure(); plot_1d(spin_system,abs(spectrum),parameters)`.

### Key state/data transformations

- Lines 11: computes `sys.magnet` using `sys.magnet=0.346`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'1H','1H','1H','E'}`.
- Lines 15: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0 0 0 2.0026}`.
- Lines 16: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(4,4)`.
- Lines 17: computes `inter.coupling.scalar{1,4}` using `inter.coupling.scalar{1,4}=-64.4e6`.
- Lines 18: computes `inter.coupling.scalar{2,4}` using `inter.coupling.scalar{2,4}=-64.4e6`.
- Lines 19: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=-64.4e6`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 24: computes `bas.sym_group` using `bas.sym_group={'S3'}`.
- Lines 25: computes `bas.sym_spins` using `bas.sym_spins={[1 2 3]}`.
- Lines 28: computes `parameters.offset` using `parameters.offset=0`.
- Lines 29: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 30: computes `parameters.sweep` using `parameters.sweep=1.2e8`.
- Lines 31: computes `parameters.tau` using `parameters.tau=100e-9`.
- Lines 32: computes `parameters.zerofill` using `parameters.zerofill=4096`.
- Lines 33: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 34: computes `parameters.axis_units` using `parameters.axis_units='MHz'`.

## Implementation structure

- Mims ENDOR spectrum of a methyl radical in liquid state.
- Magnetic parameters taken from a DFT calculation.
- Calculation time: seconds
- Magnet field
- Spin system and interactions
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Crude apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
