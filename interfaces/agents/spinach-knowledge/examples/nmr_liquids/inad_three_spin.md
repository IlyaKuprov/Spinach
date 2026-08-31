# examples/nmr_liquids/inad_three_spin.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/inad_three_spin.m`
- Signature: `inad_three_spin()`
- Total lines: 50

## Purpose

INADEQUATE spectrum of a three-spin system with J-coupling between two spins only. The sequence selects double-quantum coherence from coupled 13C pairs and converts it back for detection. Calculation time: seconds

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Magnet field; implemented by `sys.magnet=9.4`.
- Lines 15-16: Spin system and interactions; implemented by `sys.isotopes={'13C','13C','13C'}`.
- Lines 21-22: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 25-26: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 29-30: Sequence parameters; implemented by `parameters.spins={'13C'}`.
- Lines 39-40: Simulation; implemented by `fid=liquid(spin_system,@inadequate,parameters,'nmr')`.
- Lines 42-43: Processing; implemented by `fid=apodisation(spin_system,fid,{{'exp',5}})`.
- Lines 46-47: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'13C','13C','13C'}`.
- Lines 17: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={10 15 20}`.
- Lines 18: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=55`.
- Lines 19: computes `inter.coupling.scalar{3,3}` using `inter.coupling.scalar{3,3}=0`.
- Lines 22: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 23: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 26: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 30: computes `parameters.spins` using `parameters.spins={'13C'}`.
- Lines 31: computes `parameters.J` using `parameters.J=55`.
- Lines 32: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 33: computes `parameters.offset` using `parameters.offset=1800`.
- Lines 34: computes `parameters.sweep` using `parameters.sweep=5000`.
- Lines 35: computes `parameters.npoints` using `parameters.npoints=4096`.
- Lines 36: computes `parameters.zerofill` using `parameters.zerofill=16384`.
- Lines 37: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 40: computes `fid` using `fid=liquid(spin_system,@inadequate,parameters,'nmr')`.
- Lines 44: computes `spectrum` using `spectrum=fftshift(fft(fid,parameters.zerofill))`.

## Implementation structure

- INADEQUATE spectrum of a three-spin system with J-coupling between
- two spins only. The sequence selects double-quantum coherence
- from coupled 13C pairs and converts it back for detection.
- Calculation time: seconds
- Magnet field
- Spin system and interactions
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Processing
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
