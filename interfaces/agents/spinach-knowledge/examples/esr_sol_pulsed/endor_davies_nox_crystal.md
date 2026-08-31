# examples/esr_sol_pulsed/endor_davies_nox_crystal.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/endor_davies_nox_crystal.m`
- Signature: `endor_davies_nox_crystal()`
- Total lines: 116

## Purpose

Davies ENDOR simulation for a nitroxide radical at a single orientation. Soft pulses are simulated using Fokker-Planck formalism. Calculation time: minutes

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The file uses a Fokker-Planck-style enlarged state space in which spatial or orientational coordinates are promoted to extra dimensions and coupled to spin dynamics through differential operators.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 11-12: Isotopes; implemented by `sys.isotopes={'E','14N'}`.
- Lines 14-15: Magnet field; implemented by `sys.magnet=3.35`.
- Lines 17-18: Interactions; implemented by `inter.zeeman.matrix=cell(1,2)`.
- Lines 27-28: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 31-32: Relaxation theory; implemented by `inter.relaxation={'t1_t2'}`.
- Lines 38-39: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 42-45: % Stage 1: pulse-acquire ESR spectrum; implemented by `parameters.spins={'E'}`.
- Lines 44-45: Sequence parameters; implemented by `parameters.spins={'E'}`.
- Lines 58-59: Simulation; implemented by `fid=crystal(spin_system,@acquire,parameters,'esr')`.
- Lines 61-62: Apodisation; implemented by `fid=apodisation(spin_system,fid,{{'none'}})`.
- Lines 64-65: Fourier transform; implemented by `spectrum=fftshift(fft(fid,parameters.zerofill))`.
- Lines 67-68: Plotting; implemented by `kfigure(); scale_figure([1.5 1.5]); subplot(2,2,1)`.
- Lines 72-75: % Stage 2: ENDOR, left line; implemented by `clear('parameters')`.
- Lines 74-75: Sequence parameters; implemented by `clear('parameters')`.
- Lines 83-84: Electron pulse parameters; implemented by `parameters.e_rnk=2`.
- Lines 89-90: Nucleus pulse parameters; implemented by `parameters.n_rnk=3`.
- Lines 96-97: ENDOR simulation -three electron lines; implemented by `parameters.e_frq=+93e6`.
- Lines 104-105: Plotting; implemented by `subplot(2,2,2); plot(parameters.n_frq/1e6,real(answer_a))`.

### Key state/data transformations

- Lines 12: computes `sys.isotopes` using `sys.isotopes={'E','14N'}`.
- Lines 15: computes `sys.magnet` using `sys.magnet=3.35`.
- Lines 18: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,2)`.
- Lines 19: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=[2.01045 0.00000 0.00000`.
- Lines 22: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(2,2)`.
- Lines 23: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=[1.2356 0.0000 0.6322`.
- Lines 28: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 29: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 32: computes `inter.relaxation` using `inter.relaxation={'t1_t2'}`.
- Lines 33: computes `inter.rlx_keep` using `inter.rlx_keep='diagonal'`.
- Lines 34: computes `inter.equilibrium` using `inter.equilibrium='zero'`.
- Lines 35: computes `inter.r1_rates` using `inter.r1_rates={20e6 0.5e6}`.
- Lines 36: computes `inter.r2_rates` using `inter.r2_rates={20e6 0.5e6}`.
- Lines 39: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 45: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 46: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'L+','E')`.
- Lines 47: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 48: computes `parameters.decouple` using `parameters.decouple={}`.

## Implementation structure

- Davies ENDOR simulation for a nitroxide radical at a single
- orientation. Soft pulses are simulated using Fokker-Planck
- formalism.
- Calculation time: minutes
- Isotopes
- Magnet field
- Interactions
- Basis set
- Relaxation theory
- Spinach housekeeping
- % Stage 1: pulse-acquire ESR spectrum
- Sequence parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `crystal()`, `apodisation()`, `fftshift()`, `kfigure()`, `scale_figure()`, `subplot()`, `plot_1d()`, `ktitle()`, `clear()`, `kxlabel()`, `kylabel()`.
