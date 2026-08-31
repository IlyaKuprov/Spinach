# examples/esr_sol_pulsed/hard_3_pulse_echo_cu.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hard_3_pulse_echo_cu.m`
- Signature: `hard_3_pulse_echo_cu()`
- Total lines: 57

## Purpose

Three-pulse DEER echo on a Cu(II)-NO two electron system at X-band. The calculation is done by brute-force time propagation and numerical powder averaging in Liouville space. Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 12-13: Spin system parameters; implemented by `sys.magnet=0.33`.
- Lines 20-21: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 24-25: Disable trajectory level SSR algorithms; implemented by `sys.disable={'trajlevel'}`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Sequence parameters; implemented by `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 45-46: Pulse sequence; implemented by `echo=powder(spin_system,@deer_3p_hard_echo,parameters,'deer')`.
- Lines 48-50: Build the time axis; implemented by `time_axis=linspace(-parameters.tc/2, +parameters.tc/2,parameters.nsteps+1)`.
- Lines 52-53: Plotting; implemented by `kfigure(); plot(1e6*time_axis,imag(echo)); kgrid`.

### Key state/data transformations

- Lines 13: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 14: computes `sys.isotopes` using `sys.isotopes={'E','E'}`.
- Lines 15: computes `inter.zeeman.eigs` using `inter.zeeman.eigs={[2.056, 2.056, 2.205]`.
- Lines 17: computes `inter.zeeman.euler` using `inter.zeeman.euler={[0 0 0]; [0 0 0]}`.
- Lines 18: computes `inter.coordinates` using `inter.coordinates={[0 0 0]; [20 0 0]}`.
- Lines 21: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 22: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 25: computes `sys.disable` using `sys.disable={'trajlevel'}`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 33: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 34-35: computes `parameters.ex_prob` using `parameters.ex_prob=(operator(spin_system,{'L+'},{1})+ operator(spin_system,{'L-'},{1}))/2`.
- Lines 36-37: computes `parameters.ex_pump` using `parameters.ex_pump=(operator(spin_system,{'L+'},{2})+ operator(spin_system,{'L-'},{2}))/2`.
- Lines 38: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 39: computes `parameters.ta` using `parameters.ta=2.0e-7`.
- Lines 40: computes `parameters.tb` using `parameters.tb=1.0e-7`.
- Lines 41: computes `parameters.tc` using `parameters.tc=2.5e-8`.
- Lines 42: computes `parameters.nsteps` using `parameters.nsteps=256`.

## Implementation structure

- Three-pulse DEER echo on a Cu(II)-NO two electron system at X-band.
- The calculation is done by brute-force time propagation and numerical
- powder averaging in Liouville space.
- Calculation time: seconds
- Spin system parameters
- Basis set
- Disable trajectory level SSR algorithms
- Spinach housekeeping
- Sequence parameters
- Pulse sequence
- Build the time axis
- Plotting

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `kfigure()`, `kxlabel()`.
