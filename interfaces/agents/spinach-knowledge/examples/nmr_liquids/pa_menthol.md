# examples/nmr_liquids/pa_menthol.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_liquids/pa_menthol.m`
- Signature: `pa_menthol()`
- Total lines: 56

## Purpose

Menthol NMR spectrum from Damien Jeannerat, including the effect of bad Z1 and Z2 magnet shims. Calculation time: minutes.

## Physical / mathematical content

- Liquid-state NMR examples. The physics is scalar-coupling-mediated coherence transfer in weakly or moderately coupled spin systems, often in Liouville space. Typical mechanisms include INEPT-style polarisation transfer, J-refocusing, phase cycling, indirect evolution, and multidimensional detection.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 10-11: System and interaction specification; implemented by `load('menthol.mat','sys','inter')`.
- Lines 13-14: Formalism and basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 22-23: Algorithms; implemented by `sys.enable={'greedy'}`.
- Lines 25-26: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 29-30: Sequence parameters -1H; implemented by `parameters.spins={'1H'}`.
- Lines 41-42: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 44-45: Gaussian apodisation and then bad shims; implemented by `fid=apodisation(spin_system,fid,{{'gauss',5}})`.
- Lines 49-50: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 52-53: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 14: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 15: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 16: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 17: computes `bas.space_level` using `bas.space_level=1`.
- Lines 18: computes `bas.projections` using `bas.projections=+1`.
- Lines 19: computes `bas.sym_group` using `bas.sym_group={'S3','S3','S3'}`.
- Lines 20: computes `bas.sym_spins` using `bas.sym_spins={[4 5 6],[9 10 11],[12 13 14]}`.
- Lines 23: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 26: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 30: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 31: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 32: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 33: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 34: computes `parameters.offset` using `parameters.offset=1000`.
- Lines 35: computes `parameters.sweep` using `parameters.sweep=2000`.
- Lines 36: computes `parameters.npoints` using `parameters.npoints=8192`.
- Lines 37: computes `parameters.zerofill` using `parameters.zerofill=65536`.
- Lines 38: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.

## Implementation structure

- Menthol NMR spectrum from Damien Jeannerat, including the
- effect of bad Z1 and Z2 magnet shims.
- Calculation time: minutes.
- System and interaction specification
- Formalism and basis set
- Algorithms
- Spinach housekeeping
- Sequence parameters -1H
- Simulation
- Gaussian apodisation and then bad shims
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `load()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
