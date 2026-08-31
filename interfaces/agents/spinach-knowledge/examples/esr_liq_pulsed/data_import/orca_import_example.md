# examples/esr_liq_pulsed/data_import/orca_import_example.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_liq_pulsed/data_import/orca_import_example.m`
- Signature: `orca_import_example()`
- Total lines: 65

## Purpose

Methyl radical simulation, ORCA import. The uncommon signal intensity pattern comes from g-HFC cross-correlation.

## Physical / mathematical content

- Liquid-state ESR examples. The dominant physics is electron Zeeman interaction, hyperfine coupling, relaxation broadening, and pulse-acquire or ENDOR-type detection in fast tumbling systems.
- The relaxation model is Redfield-type perturbation theory: fluctuating interactions enter through correlation functions or spectral densities and generate a linear relaxation superoperator.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: System properties (vacuum DFT calculation); implemented by `props=oparse('orca_methyl_radical.out')`.
- Lines 11-12: Isotopes; implemented by `sys.isotopes={'E','1H','1H','1H'}`.
- Lines 14-15: Zeeman interactions; implemented by `inter.zeeman.matrix=cell(1,4)`.
- Lines 18-19: Hyperfine couplings; implemented by `inter.coupling.matrix=cell(4,4)`.
- Lines 24-25: Magnet induction; implemented by `sys.magnet=0.339`.
- Lines 27-28: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Relaxation theory; implemented by `inter.relaxation={'redfield'}`.
- Lines 37-38: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 41-42: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 53-54: Simulation; implemented by `fid=liquid(spin_system,@acquire,parameters,'esr')`.
- Lines 56-57: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 59-60: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 62-63: Plotting; implemented by `kfigure(); plot_1d(spin_system,real(spectrum),parameters)`.

### Key state/data transformations

- Lines 9: computes `props` using `props=oparse('orca_methyl_radical.out')`.
- Lines 12: computes `sys.isotopes` using `sys.isotopes={'E','1H','1H','1H'}`.
- Lines 15: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,4)`.
- Lines 16: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=props.g_tensor.matrix`.
- Lines 19: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(4,4)`.
- Lines 20: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=1e6*gauss2mhz(props.hfc.full.matrix{2})`.
- Lines 21: computes `inter.coupling.matrix{1,3}` using `inter.coupling.matrix{1,3}=1e6*gauss2mhz(props.hfc.full.matrix{3})`.
- Lines 22: computes `inter.coupling.matrix{1,4}` using `inter.coupling.matrix{1,4}=1e6*gauss2mhz(props.hfc.full.matrix{4})`.
- Lines 25: computes `sys.magnet` using `sys.magnet=0.339`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `inter.relaxation` using `inter.relaxation={'redfield'}`.
- Lines 33: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 34: computes `inter.rlx_keep` using `inter.rlx_keep='secular'`.
- Lines 35: computes `inter.tau_c` using `inter.tau_c={5e-10}`.
- Lines 38: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 42: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 43: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','E')`.

## Implementation structure

- Methyl radical simulation, ORCA import. The uncommon signal
- intensity pattern comes from g-HFC cross-correlation.
- System properties (vacuum DFT calculation)
- Isotopes
- Zeeman interactions
- Hyperfine couplings
- Magnet induction
- Basis set
- Relaxation theory
- Spinach housekeeping
- Sequence parameters
- Simulation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `oparse()`, `gauss2mhz()`, `create()`, `basis()`, `state()`, `liquid()`, `apodisation()`, `fftshift()`, `kfigure()`, `plot_1d()`.
