# examples/kinetics/mas_exchange_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/kinetics/mas_exchange_1.m`
- Signature: `mas_exchange_1()`
- Total lines: 69

## Purpose

Two-site position exchange for a deuterium nucleus. The sites differ in the chemical shift and the orientation of the quad- rupolar tensor. Calculation time: seconds.

## Physical / mathematical content

- Chemical-kinetics examples. The files couple spin dynamics to exchange, pumping, or nonlinear reaction networks represented by kinetic generators in Liouville space.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: System specification; implemented by `sys.magnet=9.4`.
- Lines 15-16: Spin system; implemented by `sys.isotopes={'2H','2H'}`.
- Lines 18-19: Quadrupolar interactions; implemented by `[Q1,Q2]=weblab2nqi(0.16e6,0.1,1,0,acos(1/3),2*pi/3)`.
- Lines 23-24: Chemical shifts; implemented by `inter.zeeman.scalar={0.0 3.0}`.
- Lines 26-27: Chemical exchange; implemented by `inter.chem.parts={1,2}`.
- Lines 32-33: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 36-37: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 40-41: Sequence parameters; implemented by `parameters.spins={'2H'}`.
- Lines 56-57: Acquisition; implemented by `fid=singlerot(spin_system,@acquire,parameters,'nmr')`.
- Lines 59-60: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 62-63: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 65-66: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=9.4`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'2H','2H'}`.
- Lines 19: computes `[Q1,Q2]` using `[Q1,Q2]=weblab2nqi(0.16e6,0.1,1,0,acos(1/3),2*pi/3)`.
- Lines 20: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=Q1`.
- Lines 21: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=Q2`.
- Lines 24: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0.0 3.0}`.
- Lines 27: computes `inter.chem.parts` using `inter.chem.parts={1,2}`.
- Lines 28: computes `inter.chem.rates` using `inter.chem.rates=[-2e4 2e4`.
- Lines 30: computes `inter.chem.concs` using `inter.chem.concs=[1.0 1.0]`.
- Lines 33: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 34: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 37: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 41: computes `parameters.spins` using `parameters.spins={'2H'}`.
- Lines 42: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','2H','chem')`.
- Lines 43: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','2H')`.
- Lines 44: computes `parameters.decouple` using `parameters.decouple={}`.
- Lines 45: computes `parameters.offset` using `parameters.offset=0`.
- Lines 46: computes `parameters.sweep` using `parameters.sweep=0.4e6`.

## Implementation structure

- Two-site position exchange for a deuterium nucleus. The sites
- differ in the chemical shift and the orientation of the quad-
- rupolar tensor.
- Calculation time: seconds.
- System specification
- Spin system
- Quadrupolar interactions
- Chemical shifts
- Chemical exchange
- Basis set
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `weblab2nqi()`, `acos()`, `create()`, `basis()`, `state()`, `singlerot()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
