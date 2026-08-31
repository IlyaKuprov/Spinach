# examples/esr_liq_pulsed/relaxation_parafluoronitrobenzene.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/relaxation_parafluoronitrobenzene.m`
- Signature: `relaxation_parafluoronitrobenzene()`
- Total lines: 85

## Purpose

A pulse-acquire FFT version of the EasySpin parafluoronitrobenzene test file, with acknowledgements to Stefan Stoll. The Spinach simulation is run using explicit time propagation in Liouville space, including secular Redfield relaxation superoperator. Calculation time: minutes

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Magnet field; implemented by `sys.magnet=0.33898`.
- Lines 17-18: Isotope list; implemented by `sys.isotopes={'E','14N','19F','1H','1H','1H','1H'}`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 27-28: Zeeman interactions; implemented by `inter.zeeman.eigs=cell(7,1)`.
- Lines 33-34: Spin-spin couplings; implemented by `inter.coupling.eigs=cell(7,7)`.
- Lines 49-50: Relaxation superoperator; implemented by `inter.relaxation={'redfield'}`.
- Lines 55-56: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 59-60: Experiment parameters; implemented by `parameters.spins={'E'}`.
- Lines 72-73: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'esr')`.
- Lines 75-76: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 78-79: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 81-82: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=0.33898`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'E','14N','19F','1H','1H','1H','1H'}`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 23: computes `bas.zero_quantum` using `bas.zero_quantum={'1H'}`.
- Lines 24: computes `bas.sym_group` using `bas.sym_group={'S2','S2'}`.
- Lines 25: computes `bas.sym_spins` using `bas.sym_spins={[4 5],[6 7]}`.
- Lines 28: computes `inter.zeeman.eigs` using `inter.zeeman.eigs=cell(7,1)`.
- Lines 29: computes `inter.zeeman.euler` using `inter.zeeman.euler=cell(7,1)`.
- Lines 30: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.0032 2.0012 2.0097]`.
- Lines 31: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[0 0 0]`.
- Lines 34: computes `inter.coupling.eigs` using `inter.coupling.eigs=cell(7,7)`.
- Lines 35: computes `inter.coupling.euler` using `inter.coupling.euler=cell(7,7)`.
- Lines 36: computes `inter.coupling.eigs{1,2}` using `inter.coupling.eigs{1,2}=(40.40+[24 -12 -12])*1e6`.
- Lines 37: computes `inter.coupling.eigs{1,3}` using `inter.coupling.eigs{1,3}=(22.51+[34.9 -19.8 -15])*1e6`.
- Lines 38: computes `inter.coupling.eigs{1,4}` using `inter.coupling.eigs{1,4}=[9.69 9.69 9.69]*1e6`.
- Lines 39: computes `inter.coupling.eigs{1,5}` using `inter.coupling.eigs{1,5}=[9.69 9.69 9.69]*1e6`.
- Lines 40: computes `inter.coupling.eigs{1,6}` using `inter.coupling.eigs{1,6}=[3.16 3.16 3.16]*1e6`.

## Implementation structure

- A pulse-acquire FFT version of the EasySpin parafluoronitrobenzene test
- file, with acknowledgements to Stefan Stoll.
- The Spinach simulation is run using explicit time propagation in Liouville
- space, including secular Redfield relaxation superoperator.
- Calculation time: minutes
- Magnet field
- Isotope list
- Basis set
- Zeeman interactions
- Spin-spin couplings
- Relaxation superoperator
- Spinach housekeeping

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
