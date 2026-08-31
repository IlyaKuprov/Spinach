# examples/esr_sol_pulsed/hard_3_pulse_deer_gd_2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hard_3_pulse_deer_gd_2.m`
- Signature: `hard_3_pulse_deer_gd_2()`
- Total lines: 91

## Purpose

Gadolinium(III) DEER experiment. The calculation is done by brute- force time propagation and powder averaging. Outermost ZFS transi- tion is excited by the probe pulse and the central transition is excited by the pump pulse. The pulses are assumed to be ideal. Note: gadolinium spin echo is very sharp and difficult to catch in simulations because they do not include zero-field splitting distributions found in experim

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- Signal processing is central here: the code moves between time and frequency domains, typically using FFT conventions, apodisation, zero filling, or heterodyne frequency shifts.

## Numerical / algorithmic content

- The output is processed in the Fourier domain, implying standard NMR/ESR signal-processing considerations such as acquisition bandwidth, zero filling, phase, and apodisation.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 19-20: Spin system properties; implemented by `sys.magnet=3.5`.
- Lines 29-30: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 33-34: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 37-38: Probe pulse operator -bottom transition; implemented by `sigma=pauli(2); sigma.p=[zeros(6,8); [zeros(2,6) sigma.p]]`.
- Lines 41-42: Pump pulse operator -central transition; implemented by `sigma=pauli(2); sigma.p=[zeros(3,8); [zeros(2,3) sigma.p zeros(2,3)]; zeros(3,8)]`.
- Lines 45-46: Sequence parameters; implemented by `parameters.rho0=state(spin_system,'Lz','E8')`.
- Lines 60-61: Pulse sequence; implemented by `deer=powder(spin_system,@deer_3p_hard_deer,parameters,'deer')`.
- Lines 63-64: Apodisation; implemented by `deer.hard_pulse_fid=apodisation(spin_system,deer.hard_pulse_fid,{{'exp',6}})`.
- Lines 68-69: Fourier transforms; implemented by `hard_pulse_spec=imag(fftshift(fft(deer.hard_pulse_fid,4*parameters.spectrum_nsteps)))`.
- Lines 74-75: Plotting; implemented by `kfigure(); scale_figure([1.0 2.0])`.

### Key state/data transformations

- Lines 20: computes `sys.magnet` using `sys.magnet=3.5`.
- Lines 21: computes `sys.isotopes` using `sys.isotopes={'E8','E8'}`.
- Lines 22: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.002319 2.002319}`.
- Lines 23: computes `inter.coupling.eigs{1,1}` using `inter.coupling.eigs{1,1}=[1e9 1e9 -2e9]`.
- Lines 24: computes `inter.coupling.euler{1,1}` using `inter.coupling.euler{1,1}=[0 pi/9 0]`.
- Lines 25: computes `inter.coupling.eigs{2,2}` using `inter.coupling.eigs{2,2}=[1e9 1e9 -2e9]`.
- Lines 26: computes `inter.coupling.euler{2,2}` using `inter.coupling.euler{2,2}=[0 4*pi/9 0]`.
- Lines 27: computes `inter.coordinates` using `inter.coordinates={[0 0 0]; 30*[sind(20) 0 cosd(20)]}`.
- Lines 30: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 31: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 34: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 38: computes `sigma` using `sigma=pauli(2); sigma.p=[zeros(6,8); [zeros(2,6) sigma.p]]`.
- Lines 39: computes `Ep_prob` using `Ep_prob=kron(sigma.p,speye(size(sigma.p)))`.
- Lines 43: computes `Ep_pump` using `Ep_pump=kron(speye(size(sigma.p)),sigma.p)`.
- Lines 46: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E8')`.
- Lines 47: computes `parameters.ex_prob` using `parameters.ex_prob=(Ep_prob+Ep_prob')/2`.
- Lines 48: computes `parameters.ex_pump` using `parameters.ex_pump=(Ep_pump+Ep_pump')/2`.
- Lines 49: computes `parameters.coil_prob` using `parameters.coil_prob=state(spin_system,{'L+'},{1})`.

## Implementation structure

- Gadolinium(III) DEER experiment. The calculation is done by brute-
- force time propagation and powder averaging. Outermost ZFS transi-
- tion is excited by the probe pulse and the central transition is
- excited by the pump pulse. The pulses are assumed to be ideal.
- Note: gadolinium spin echo is very sharp and difficult to catch in
- simulations because they do not include zero-field splitting
- distributions found in experimental systems.
- Calculation time: minutes.
- Spin system properties
- Basis set
- Spinach housekeeping
- Probe pulse operator -bottom transition

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `sind()`, `cosd()`, `create()`, `basis()`, `pauli()`, `speye()`, `state()`, `operator()`, `powder()`, `apodisation()`, `fftshift()`, `ft_axis()`, `kfigure()`, `scale_figure()`, `subplot()`, `kxlabel()`.
