# examples/parahydrogen/ortho_deuterium.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/parahydrogen/ortho_deuterium.m`
- Signature: `ortho_deuterium()`
- Total lines: 68

## Purpose

Ortho-deuteration simulation for acrylonitrile in Figure 1 of the paper by Natterer, Greve, and Bargon: Simulation time: seconds

## Physical / mathematical content

- Parahydrogen examples. The physical motif is highly non-Boltzmann singlet order imported from para-H2 and converted into observable nuclear magnetisation through hydrogenation, exchange, or catalytic transfer processes.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Bargon's magnet; implemented by `sys.magnet=4.697`.
- Lines 16-17: Deuteration product; implemented by `sys.isotopes={'1H','1H','2H','1H','2H'}`.
- Lines 28-29: Hilbert space; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Continuous deuteration; implemented by `options.dephasing=1`.
- Lines 39-40: Singlet and quintet on deuterium; implemented by `[S,~,Q]=deut_pair(spin_system,3,5,options)`.
- Lines 42-43: Experiment parameters; implemented by `parameters.spins={'2H'}`.
- Lines 55-56: Simulation; implemented by `fid=liquid(spin_system,@hp_acquire,parameters,'nmr')`.
- Lines 58-59: Apodisation and sign flip; implemented by `fid=apodisation(spin_system,fid,{{'exp',6}})`.
- Lines 61-62: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 64-65: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=4.697`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'1H','1H','2H','1H','2H'}`.
- Lines 18: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={1.2 1.2 1.2 2.3 2.3}`.
- Lines 19: computes `inter.coupling.scalar` using `inter.coupling.scalar=cell(5,5)`.
- Lines 20: computes `inter.coupling.scalar{1,3}` using `inter.coupling.scalar{1,3}=2.0`.
- Lines 21: computes `inter.coupling.scalar{2,3}` using `inter.coupling.scalar{2,3}=2.0`.
- Lines 22: computes `inter.coupling.scalar{4,5}` using `inter.coupling.scalar{4,5}=2.0`.
- Lines 23: computes `inter.coupling.scalar{3,4}` using `inter.coupling.scalar{3,4}=1.2`.
- Lines 24: computes `inter.coupling.scalar{1,5}` using `inter.coupling.scalar{1,5}=1.2`.
- Lines 25: computes `inter.coupling.scalar{2,5}` using `inter.coupling.scalar{2,5}=1.2`.
- Lines 26: computes `inter.coupling.scalar{3,5}` using `inter.coupling.scalar{3,5}=0.2`.
- Lines 29: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `options.dephasing` using `options.dephasing=1`.
- Lines 40: computes `[S,~,Q]` using `[S,~,Q]=deut_pair(spin_system,3,5,options)`.
- Lines 43: computes `parameters.spins` using `parameters.spins={'2H'}`.
- Lines 44: computes `parameters.rho0` using `parameters.rho0=S+Q{1}+Q{2}+Q{3}+Q{4}+Q{5}`.

## Implementation structure

- Ortho-deuteration simulation for acrylonitrile in Figure 1 of
- the paper by Natterer, Greve, and Bargon:
- Simulation time: seconds
- Bargon's magnet
- Deuteration product
- Hilbert space
- Spinach housekeeping
- Continuous deuteration
- Singlet and quintet on deuterium
- Experiment parameters
- Simulation
- Apodisation and sign flip

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `deut_pair()`, `state()`, `operator()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
