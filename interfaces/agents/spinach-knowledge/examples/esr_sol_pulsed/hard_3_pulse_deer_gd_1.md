# examples/esr_sol_pulsed/hard_3_pulse_deer_gd_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hard_3_pulse_deer_gd_1.m`
- Signature: `hard_3_pulse_deer_gd_1()`
- Total lines: 94

## Purpose

Gadolinium(III) DEER experiment at W-band using ideal pulses. Set to reproduce Figure 2b from the paper by Otting and co-authors: The calculation is done by brute-force time propagation and grid pow- der averaging. Central transitions are used on both gadolinium ions. Note: gadolinium spin echo is very sharp and difficult to catch in simulations because they do not include zero-field splitting distributions found in 

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 26-27: Spin system properties; implemented by `sys.magnet=3.5`.
- Lines 40-41: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 44-45: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 48-49: Sequence parameters; implemented by `parameters.rho0=state(spin_system,'Lz','E8')`.
- Lines 63-64: Pulse sequence; implemented by `deer=powder(spin_system,@deer_3p_hard_deer,parameters,'deer-zz')`.
- Lines 66-67: Apodisation; implemented by `deer.hard_pulse_fid=apodisation(spin_system,deer.hard_pulse_fid,{{'exp',6}})`.
- Lines 71-72: Fourier transforms; implemented by `hard_pulse_spec=imag(fftshift(fft(deer.hard_pulse_fid,4*parameters.spectrum_nsteps)))`.
- Lines 77-78: Plotting; implemented by `kfigure(); scale_figure([1.0 2.0])`.

### Key state/data transformations

- Lines 27: computes `sys.magnet` using `sys.magnet=3.5`.
- Lines 28: computes `sys.isotopes` using `sys.isotopes={'E8','E8'}`.
- Lines 29: computes `inter.zeeman.scalar{1}` using `inter.zeeman.scalar{1}=2.002319`.
- Lines 30: computes `inter.zeeman.scalar{2}` using `inter.zeeman.scalar{2}=2.002319`.
- Lines 31: computes `inter.coordinates` using `inter.coordinates={[ 0.00 0.00 0.00]`.
- Lines 33: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=[1e8 0 0`.
- Lines 36: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=[1e8 0 0`.
- Lines 41: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 42: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 45: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 49: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E8')`.
- Lines 50: computes `parameters.ex_prob` using `parameters.ex_prob=operator(spin_system,'CTx',1)`.
- Lines 51: computes `parameters.ex_pump` using `parameters.ex_pump=operator(spin_system,'CTx',2)`.
- Lines 52: computes `parameters.coil_prob` using `parameters.coil_prob=state(spin_system,{'L+'},{1})`.
- Lines 53: computes `parameters.coil_pump` using `parameters.coil_pump=state(spin_system,{'L+'},{2})`.
- Lines 54: computes `parameters.spectrum_sweep` using `parameters.spectrum_sweep=1e10`.
- Lines 55: computes `parameters.spectrum_nsteps` using `parameters.spectrum_nsteps=1024`.
- Lines 56: computes `parameters.ex_hard` using `parameters.ex_hard=operator(spin_system,'Lx','electrons')`.

## Implementation structure

- Gadolinium(III) DEER experiment at W-band using ideal pulses. Set to
- reproduce Figure 2b from the paper by Otting and co-authors:
- The calculation is done by brute-force time propagation and grid pow-
- der averaging. Central transitions are used on both gadolinium ions.
- Note: gadolinium spin echo is very sharp and difficult to catch in
- simulations because they do not include zero-field splitting
- distributions found in experimental systems.
- Note: flip-flop terms in the inter-electron dipolar interaction are
- switched off ('deer-zz') to mimic the effects of slightly dif-
- ferent pulse frequencies in the experiment.
- Calculation time: minutes.
- Spin system properties

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `apodisation()`, `fftshift()`, `ft_axis()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`, `ktitle()`.
