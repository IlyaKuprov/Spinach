# examples/esr_sol_pulsed/hard_3_pulse_echo_no.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hard_3_pulse_echo_no.m`
- Signature: `hard_3_pulse_echo_no()`
- Total lines: 58

## Purpose

DEER spin echo for a pair of nitroxide radicals at X-band. Two nit- roxide radicals are positioned at a distance of 25 Angstroms. The calculation is done by brute-force time propagation and numerical powder averaging in Liouville space. Nitroxide g-tensor data comes from http://dx.doi.org/10.1063/1.1697233 Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 14-15: Spin system properties; implemented by `sys.magnet=0.33`.
- Lines 24-25: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 28-29: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 32-33: Sequence parameters; implemented by `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 46-47: Pulse sequence; implemented by `echo=powder(spin_system,@deer_3p_hard_echo,parameters,'deer')`.
- Lines 49-51: Build the time axis; implemented by `time_axis=linspace(-parameters.tc/2, +parameters.tc/2,parameters.nsteps+1)`.
- Lines 53-54: Plotting; implemented by `kfigure(); plot(1e6*time_axis,imag(echo)); kgrid`.

### Key state/data transformations

- Lines 15: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 16: computes `sys.isotopes` using `sys.isotopes={'E','E'}`.
- Lines 17: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.0089 2.0061 2.0027]`.
- Lines 18: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[2.0089 2.0061 2.0027]`.
- Lines 19: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[1.0 2.0 3.0]`.
- Lines 20: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[3.0 1.0 2.0]`.
- Lines 21: computes `inter.coordinates` using `inter.coordinates={[ 0.00 0.00 0.00]`.
- Lines 25: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 26: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 29: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 33: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 34: computes `parameters.coil` using `parameters.coil=state(spin_system,'L+','E')`.
- Lines 35-36: computes `parameters.ex_prob` using `parameters.ex_prob=(operator(spin_system,{'L+'},{1})+ operator(spin_system,{'L-'},{1}))/2`.
- Lines 37-38: computes `parameters.ex_pump` using `parameters.ex_pump=(operator(spin_system,{'L+'},{2})+ operator(spin_system,{'L-'},{2}))/2`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 40: computes `parameters.ta` using `parameters.ta=1.0e-6`.
- Lines 41: computes `parameters.tb` using `parameters.tb=0.5e-6`.
- Lines 42: computes `parameters.tc` using `parameters.tc=0.5e-6`.

## Implementation structure

- DEER spin echo for a pair of nitroxide radicals at X-band. Two nit-
- roxide radicals are positioned at a distance of 25 Angstroms.
- The calculation is done by brute-force time propagation and numerical
- powder averaging in Liouville space. Nitroxide g-tensor data comes
- from http://dx.doi.org/10.1063/1.1697233
- Calculation time: seconds
- Spin system properties
- Basis set
- Spinach housekeeping
- Sequence parameters
- Pulse sequence
- Build the time axis

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `kfigure()`, `kxlabel()`.
