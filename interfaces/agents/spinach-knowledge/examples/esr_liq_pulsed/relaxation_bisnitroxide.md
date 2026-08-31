# examples/esr_liq_pulsed/relaxation_bisnitroxide.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/relaxation_bisnitroxide.m`
- Signature: `relaxation_bisnitroxide()`
- Total lines: 79

## Purpose

X-band pulse-acquire FFT ESR spectrum of a bisnitroxide radical, using explicit time domain simulation with Redfield relaxation supeoperator. Parameters from https://doi.org/10.1039/C8CP06819D Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Spin system; implemented by `sys.isotopes={'E','E','14N','14N'}`.
- Lines 14-15: Zeeman interactions; implemented by `inter.zeeman.eigs{1}=[2.00925 2.00605 2.00205]`.
- Lines 24-25: Couplings; implemented by `inter.coupling.eigs=cell(4,4)`.
- Lines 36-37: Magnet induction; implemented by `sys.magnet=0.35`.
- Lines 39-40: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 43-44: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 49-50: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 53-54: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 66-67: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'esr')`.
- Lines 69-70: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 72-73: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 75-76: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'E','E','14N','14N'}`.
- Lines 15: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.00925 2.00605 2.00205]`.
- Lines 16: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[0 0 0]`.
- Lines 17: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=inter.zeeman.eigs{1}`.
- Lines 18: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[123.1 129.8 -46.6]*pi/180`.
- Lines 19: computes `inter.zeeman.eigs{3}` using `inter.zeeman.eigs{3}=[0 0 0]`.
- Lines 20: computes `inter.zeeman.euler{3}` using `inter.zeeman.euler{3}=[0 0 0]`.
- Lines 21: computes `inter.zeeman.eigs{4}` using `inter.zeeman.eigs{4}=[0 0 0]`.
- Lines 22: computes `inter.zeeman.euler{4}` using `inter.zeeman.euler{4}=[0 0 0]`.
- Lines 25: computes `inter.coupling.eigs` using `inter.coupling.eigs=cell(4,4)`.
- Lines 26: computes `inter.coupling.euler` using `inter.coupling.euler=cell(4,4)`.
- Lines 27: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(4,4)`.
- Lines 28: computes `inter.coupling.eigs{1,2}` using `inter.coupling.eigs{1,2}=[17.5,17.5,-35]*1e6`.
- Lines 29: computes `inter.coupling.eigs{1,3}` using `inter.coupling.eigs{1,3}=[18,17,103]*1e6`.
- Lines 30: computes `inter.coupling.eigs{2,4}` using `inter.coupling.eigs{2,4}=inter.coupling.eigs{1,3}`.
- Lines 31: computes `inter.coupling.euler{1,2}` using `inter.coupling.euler{1,2}=[-174 74 0]*pi/180`.
- Lines 32: computes `inter.coupling.euler{1,3}` using `inter.coupling.euler{1,3}=[0 0 0]*pi/180`.
- Lines 33: computes `inter.coupling.euler{2,4}` using `inter.coupling.euler{2,4}=inter.zeeman.euler{2}`.

## Implementation structure

- X-band pulse-acquire FFT ESR spectrum of a bisnitroxide radical,
- using explicit time domain simulation with Redfield relaxation
- supeoperator. Parameters from https://doi.org/10.1039/C8CP06819D
- Calculation time: seconds
- Spin system
- Zeeman interactions
- Couplings
- Magnet induction
- Basis set
- Relaxation theory
- Spinach housekeeping
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
