# examples/esr_liq_pulsed/relaxation_fremysalt.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/relaxation_fremysalt.m`
- Signature: `relaxation_fremysalt()`
- Total lines: 72

## Purpose

Pulse-acquire FFT ESR version of the EasySpin Fremy salt test file, with acknowledgements to Stefan Stoll. The Spinach simulation is run using explicit time propagation in Liouville space with Redfield relaxation superoperator. Set to reproduce Figure 3a from Calculation time: seconds

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: General layout; implemented by `sys.magnet=0.33`.
- Lines 22-23: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 26-27: Interactions; implemented by `inter.zeeman.eigs=cell(2,1)`.
- Lines 36-37: Relaxation superoperator; implemented by `inter.relaxation={'redfield'}`.
- Lines 42-43: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 46-47: Experiment parameters; implemented by `parameters.spins={'E'}`.
- Lines 59-60: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'esr')`.
- Lines 62-63: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 65-66: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 68-69: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 19: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 20: computes `sys.isotopes` using `sys.isotopes={'E','14N'}`.
- Lines 23: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 24: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 27: computes `inter.zeeman.eigs` using `inter.zeeman.eigs=cell(2,1)`.
- Lines 28: computes `inter.zeeman.euler` using `inter.zeeman.euler=cell(2,1)`.
- Lines 29: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.00785 2.00590 2.00265]`.
- Lines 30: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[0 0 0]`.
- Lines 31: computes `inter.coupling.eigs` using `inter.coupling.eigs=cell(2,2)`.
- Lines 32: computes `inter.coupling.euler` using `inter.coupling.euler=cell(2,2)`.
- Lines 33: computes `inter.coupling.eigs{1,2}` using `inter.coupling.eigs{1,2}=[15.4137 14.0125 80.4316]*1e6`.
- Lines 34: computes `inter.coupling.euler{1,2}` using `inter.coupling.euler{1,2}=[0 0 0]`.
- Lines 37: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 38: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 39: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 40: computes `inter.tau_c` using `inter.tau_c={8e-10}`.
- Lines 43: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 47: computes `parameters.spins` using `parameters.spins={'E'}`.

## Implementation structure

- Pulse-acquire FFT ESR version of the EasySpin Fremy salt test
- file, with acknowledgements to Stefan Stoll.
- The Spinach simulation is run using explicit time propagation
- in Liouville space with Redfield relaxation superoperator.
- Set to reproduce Figure 3a from
- Calculation time: seconds
- General layout
- Basis set
- Interactions
- Relaxation superoperator
- Spinach housekeeping
- Experiment parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
