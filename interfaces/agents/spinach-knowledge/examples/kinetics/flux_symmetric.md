# examples/kinetics/flux_symmetric.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/flux_symmetric.m`
- Signature: `flux_symmetric()`
- Total lines: 52

## Purpose

Two-spin symmetric magnetization flux problem. Calculation time: seconds.

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: System specification; implemented by `sys.magnet=14.1`.
- Lines 18-19: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 22-23: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 26-27: Sequence parameters; implemented by `parameters.spins={'1H'}`.
- Lines 39-40: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 42-43: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 45-46: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 48-49: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 10: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 11: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 12: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0 3.0}`.
- Lines 13: computes `inter.chem.flux_rate` using `inter.chem.flux_rate=zeros(2)`.
- Lines 14: computes `inter.chem.flux_rate(1,2)` using `inter.chem.flux_rate(1,2)=2e3`.
- Lines 15: computes `inter.chem.flux_rate(2,1)` using `inter.chem.flux_rate(2,1)=2e3`.
- Lines 16: computes `inter.chem.flux_type` using `inter.chem.flux_type='intermolecular'`.
- Lines 19: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 20: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 23: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 27: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 28: computes `parameters.rho0` using `parameters.rho0=1.0*state(spin_system,{'L+'},{1})+`.
- Lines 30: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H')`.
- Lines 31: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 32: computes `parameters.offset` using `parameters.offset=900`.
- Lines 33: computes `parameters.sweep` using `parameters.sweep=5000`.
- Lines 34: computes `parameters.npoints` using `parameters.npoints=512`.
- Lines 35: computes `parameters.zerofill` using `parameters.zerofill=1024`.

## Implementation structure

- Two-spin symmetric magnetization flux problem.
- Calculation time: seconds.
- System specification
- Basis set
- Spinach housekeeping
- Sequence parameters
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
