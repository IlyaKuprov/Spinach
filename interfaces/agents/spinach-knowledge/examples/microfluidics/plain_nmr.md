# examples/microfluidics/plain_nmr.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/microfluidics/plain_nmr.m`
- Signature: `plain_nmr()`
- Total lines: 52

## Purpose

NMR spectrum of the reaction mixture in the absence of chemical kinetics and spatial dynamics.

## Physical / mathematical content

- Microfluidics examples. The coupled model is spin dynamics plus advection-diffusion-reaction transport on a mesh or regular grid. Numerical issues include finite-difference operators, mesh interpolation, and coupled reaction-flow evolution.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Import Diels-Alder cycloaddition; implemented by `[sys,inter,bas]=dac_reaction()`.
- Lines 14-15: Equal concentrations, no solvent; implemented by `inter.chem.concs=[1 1 1 1 0]`.
- Lines 17-18: Magnet field; implemented by `sys.magnet=14.1`.
- Lines 20-21: Greedy parallelisation; implemented by `sys.enable={'greedy'}`.
- Lines 23-24: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 27-28: Sequence parameters -1H; implemented by `parameters.spins={'1H'}`.
- Lines 39-40: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 42-43: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 45-46: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 48-49: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 12: computes `[sys,inter,bas]` using `[sys,inter,bas]=dac_reaction()`.
- Lines 15: computes `inter.chem.concs` using `inter.chem.concs=[1 1 1 1 0]`.
- Lines 18: computes `sys.magnet` using `sys.magnet=14.1`.
- Lines 21: computes `sys.enable` using `sys.enable={'greedy'}`.
- Lines 24: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 28: computes `parameters.spins` using `parameters.spins={'1H'}`.
- Lines 29: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','1H','chem')`.
- Lines 30: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','1H','chem')`.
- Lines 31: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 32: computes `parameters.offset` using `parameters.offset=2328`.
- Lines 33: computes `parameters.sweep` using `parameters.sweep=3500`.
- Lines 34: computes `parameters.npoints` using `parameters.npoints=4096`.
- Lines 35: computes `parameters.zerofill` using `parameters.zerofill=16384`.
- Lines 36: computes `parameters.axis_units` using `parameters.axis_units='ppm'`.
- Lines 37: computes `parameters.invert_axis` using `parameters.invert_axis=1`.
- Lines 40: computes `fid` using `fid=liquid(spin_system,@acquire,parameters,'nmr')`.
- Lines 46: computes `spectrum` using `spectrum=fftshift(fft(fid,parameters.zerofill))`.

## Implementation structure

- NMR spectrum of the reaction mixture in the absence of
- chemical kinetics and spatial dynamics.
- Import Diels-Alder cycloaddition
- Equal concentrations, no solvent
- Magnet field
- Greedy parallelisation
- Spinach housekeeping
- Sequence parameters -1H
- Simulation
- Apodisation
- Fourier transform
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `dac_reaction()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
