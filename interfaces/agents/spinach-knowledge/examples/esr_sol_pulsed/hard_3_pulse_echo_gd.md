# examples/esr_sol_pulsed/hard_3_pulse_echo_gd.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hard_3_pulse_echo_gd.m`
- Signature: `hard_3_pulse_echo_gd()`
- Total lines: 68

## Purpose

Gadolinium(III) DEER echo experiment. The calculation is done by brute-force time propagation and powder averaging. Outermost ZFS transition is excited by the probe pulse and the central transi- tion is excited by the pump pulse. Pulses are assumed to be hard. Note: gadolinium spin echo is very sharp and difficult to catch in simulations because they do not include zero-field splitting distributions found in experime

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 18-19: Spin system properties; implemented by `sys.magnet=3.5; D=0.56e9`.
- Lines 28-29: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 32-33: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 36-37: Probe pulse operator; implemented by `sigma=pauli(2); sigma.p=[zeros(6,8); [zeros(2,6) sigma.p]]`.
- Lines 40-41: Pump pulse operator; implemented by `sigma=pauli(2); sigma.p=[zeros(3,8); [zeros(2,3) sigma.p zeros(2,3)]; zeros(3,8)]`.
- Lines 44-45: Sequence parameters; implemented by `parameters.rho0=state(spin_system,'Lz','E8')`.
- Lines 56-57: Pulse sequence; implemented by `echo=powder(spin_system,@deer_3p_hard_echo,parameters,'deer')`.
- Lines 59-61: Build the time axis; implemented by `time_axis=linspace(-parameters.tc/2, +parameters.tc/2,parameters.nsteps+1)`.
- Lines 63-64: Plotting; implemented by `kfigure(); plot(1e6*time_axis,imag(echo)); kgrid`.

### Key state/data transformations

- Lines 19: computes `sys.magnet` using `sys.magnet=3.5; D=0.56e9`.
- Lines 20: computes `sys.isotopes` using `sys.isotopes={'E8','E8'}`.
- Lines 21: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={2.002319 2.002319}`.
- Lines 22: computes `inter.coordinates` using `inter.coordinates={[0 0 0]; 29.5*[0 1 0]}`.
- Lines 23: computes `inter.coupling.eigs{1,1}` using `inter.coupling.eigs{1,1}=[D D -2*D]/3`.
- Lines 24: computes `inter.coupling.eigs{2,2}` using `inter.coupling.eigs{2,2}=[D D -2*D]/3`.
- Lines 25: computes `inter.coupling.euler{1,1}` using `inter.coupling.euler{1,1}=[0 0 0]`.
- Lines 26: computes `inter.coupling.euler{2,2}` using `inter.coupling.euler{2,2}=[0 pi/2 0]`.
- Lines 29: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 30: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 33: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 37: computes `sigma` using `sigma=pauli(2); sigma.p=[zeros(6,8); [zeros(2,6) sigma.p]]`.
- Lines 38: computes `Ep_prob` using `Ep_prob=kron(sigma.p,speye(size(sigma.p)))`.
- Lines 42: computes `Ep_pump` using `Ep_pump=kron(speye(size(sigma.p)),sigma.p)`.
- Lines 45: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E8')`.
- Lines 46: computes `parameters.ex_prob` using `parameters.ex_prob=(Ep_prob+Ep_prob')/2`.
- Lines 47: computes `parameters.ex_pump` using `parameters.ex_pump=(Ep_pump+Ep_pump')/2`.
- Lines 48: computes `parameters.coil` using `parameters.coil=state(spin_system,{'L+'},{1})`.

## Implementation structure

- Gadolinium(III) DEER echo experiment. The calculation is done by
- brute-force time propagation and powder averaging. Outermost ZFS
- transition is excited by the probe pulse and the central transi-
- tion is excited by the pump pulse. Pulses are assumed to be hard.
- Note: gadolinium spin echo is very sharp and difficult to catch in
- simulations because they do not include zero-field splitting
- distributions found in experimental systems.
- Calculation time: seconds
- Spin system properties
- Basis set
- Spinach housekeeping
- Probe pulse operator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `pauli()`, `speye()`, `state()`, `powder()`, `kfigure()`, `kxlabel()`.
