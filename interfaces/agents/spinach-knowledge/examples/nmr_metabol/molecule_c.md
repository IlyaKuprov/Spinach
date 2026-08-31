# examples/nmr_metabol/molecule_c.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/nmr_metabol/molecule_c.m`
- Signature: `molecule_c()`
- Total lines: 47

## Purpose

1H NMR spectrum of a molecule from the GISSMO database. Calculation time: seconds

## Physical / mathematical content

- Metabolomics NMR examples. These files apply liquid-state NMR simulation workflows to small-molecule mixtures, spectral assignment, concentration inference, and database-style metabolite spin-system definitions.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Import GISSMO dataset; implemented by `[sys,inter]=gissmo2spinach('molecule_c.xml',1)`.
- Lines 12-13: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 18-19: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 22-23: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 34-35: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 37-38: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'gauss',10}})`.
- Lines 40-41: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 43-44: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 10: computes `[sys,inter]` using `[sys,inter]=gissmo2spinach('molecule_c.xml',1)`.
- Lines 13: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 14: computes `bas.approximation` using `bas.approximation='IK-2'`.
- Lines 15: computes `bas.connectivity` using `bas.connectivity='scalar_couplings'`.
- Lines 16: computes `bas.space_level` using `bas.space_level=1`.
- Lines 19: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 23: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 24: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H')`.
- Lines 25: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 26: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 27: computes `parameters.offset` using `parameters.offset=2500`.
- Lines 28: computes `parameters.sweep` using `parameters.sweep=5000`.
- Lines 29: computes `parameters.npoints` using `parameters.npoints=4096`.
- Lines 30: computes `parameters.zerofill` using `parameters.zerofill=16536`.
- Lines 31: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 32: computes `parameters.invert_axis` using `parameters.invert_axis=1`.
- Lines 35: computes `fid` using `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 41: computes `spectrum` using `spectrum=fftshift(fft(fid,parameters.zerofill))`.

## Implementation structure

- 1H NMR spectrum of a molecule from the GISSMO database.
- Calculation time: seconds
- Import GISSMO dataset
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gissmo2spinach()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
