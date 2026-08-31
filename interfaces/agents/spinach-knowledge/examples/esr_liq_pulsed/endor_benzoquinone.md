# examples/esr_liq_pulsed/endor_benzoquinone.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/endor_benzoquinone.m`
- Signature: `endor_benzoquinone()`
- Total lines: 57

## Purpose

CW ENDOR on 2-methoxy-1,4-benzoquinone radical in liquid state. Set to reproduce Figure 2 in http://dx.doi.org/10.1002/mrc.1260280313 Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Magnet field; implemented by `sys.magnet=0.33`.
- Lines 14-15: Isotopes and interactions; implemented by `sys.isotopes={'E','1H','1H','1H','1H','1H','1H'}`.
- Lines 25-26: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Sequence parameters; implemented by `parameters.offset=0`.
- Lines 40-41: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 44-45: Simulation; implemented by `fid=liquid(spin_system,@endor_cw,parameters,'esr')`.
- Lines 47-48: Crude apodization; implemented by `fid=apodisation(spin_system,fid-mean(fid),{{'kaiser',20}})`.
- Lines 50-51: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 53-54: Plotting; implemented by `kfigure(); plot_1d(spin_system,-abs(spectrum),parameters)`.

### Key state/data transformations

- Lines 12: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 15: computes `sys.isotopes` using `sys.isotopes={'E','1H','1H','1H','1H','1H','1H'}`.
- Lines 16: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.004577 0 0 0 0 0 0}`.
- Lines 17: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(7,7)`.
- Lines 18: computes `inter.coupling.scalar{2,1}` using `inter.coupling.scalar{2,1}=mt2hz(0.08)`.
- Lines 19: computes `inter.coupling.scalar{3,1}` using `inter.coupling.scalar{3,1}=mt2hz(0.08)`.
- Lines 20: computes `inter.coupling.scalar{4,1}` using `inter.coupling.scalar{4,1}=mt2hz(0.08)`.
- Lines 21: computes `inter.coupling.scalar{5,1}` using `inter.coupling.scalar{5,1}=mt2hz(-0.059)`.
- Lines 22: computes `inter.coupling.scalar{6,1}` using `inter.coupling.scalar{6,1}=mt2hz(-0.364)`.
- Lines 23: computes `inter.coupling.scalar{7,1}` using `inter.coupling.scalar{7,1}=mt2hz(-0.204)`.
- Lines 26: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 27: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 28: computes `bas.sym_group` using `bas.sym_group={'S3'}`.
- Lines 29: computes `bas.sym_spins` using `bas.sym_spins={[2 3 4]}`.
- Lines 32: computes `parameters.offset` using `parameters.offset=0`.
- Lines 33: computes `parameters.sweep` using `parameters.sweep=50e6`.
- Lines 34: computes `parameters.npoints` using `parameters.npoints=1024`.
- Lines 35: computes `parameters.zerofill` using `parameters.zerofill=4096`.

## Implementation structure

- CW ENDOR on 2-methoxy-1,4-benzoquinone radical in liquid state. Set to
- reproduce Figure 2 in http://dx.doi.org/10.1002/mrc.1260280313
- Calculation time: seconds
- Magnet field
- Isotopes and interactions
- Basis set
- Sequence parameters
- Spinach housekeeping
- Simulation
- Crude apodization
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mt2hz()`, `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
