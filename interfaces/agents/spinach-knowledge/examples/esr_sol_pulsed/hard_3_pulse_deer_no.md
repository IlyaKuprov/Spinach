# examples/esr_sol_pulsed/hard_3_pulse_deer_no.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hard_3_pulse_deer_no.m`
- Signature: `hard_3_pulse_deer_no()`
- Total lines: 56

## Purpose

Nitroxide spin label DEER experiment at X-band. Two nitroxide radicals are positioned at a distance of 25 Angstroms. The numerical calculation is done by brute-force time propaga- tion and numerical powder averaging, including g-factor orien- tation effects on the dipolar coupling. Nitroxide g-tensor is from http://dx.doi.org/10.1063/1.1697233 Calculation time: seconds

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Spin system properties; implemented by `sys.magnet=0.33`.
- Lines 26-27: Basis set; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 30-31: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 34-35: Sequence parameters; implemented by `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 45-46: Pulse sequence; implemented by `deer=powder(spin_system,@deer_3p_hard_deer,parameters,'deer')`.
- Lines 48-49: Build the time axis; implemented by `time_axis=linspace(0,parameters.stepsize*parameters.nsteps,parameters.nsteps+1)`.
- Lines 51-52: Plotting; implemented by `kfigure(); plot(1e6*time_axis,imag(deer.deer_trace))`.

### Key state/data transformations

- Lines 17: computes `sys.magnet` using `sys.magnet=0.33`.
- Lines 18: computes `sys.isotopes` using `sys.isotopes={'E','E'}`.
- Lines 19: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.0089 2.0061 2.0027]`.
- Lines 20: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[2.0089 2.0061 2.0027]`.
- Lines 21: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[1.0 2.0 3.0]`.
- Lines 22: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[3.0 1.0 2.0]`.
- Lines 23: computes `inter.coordinates` using `inter.coordinates={[ 0.00 0.00 0.00]`.
- Lines 27: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 28: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 31: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 35: computes `parameters.rho0` using `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 36: computes `parameters.coil_prob` using `parameters.coil_prob=state(spin_system,{'L-'},{1})`.
- Lines 37: computes `parameters.stepsize` using `parameters.stepsize=1e-8`.
- Lines 38: computes `parameters.nsteps` using `parameters.nsteps=100`.
- Lines 39: computes `parameters.spins` using `parameters.spins={'E'}`.
- Lines 40: computes `parameters.ex_prob` using `parameters.ex_prob=operator(spin_system,{'Lx'},{1})`.
- Lines 41: computes `parameters.ex_pump` using `parameters.ex_pump=operator(spin_system,{'Lx'},{2})`.
- Lines 42: computes `parameters.output` using `parameters.output='brief'`.

## Implementation structure

- Nitroxide spin label DEER experiment at X-band. Two nitroxide radicals
- are positioned at a distance of 25 Angstroms.
- The numerical calculation is done by brute-force time propaga-
- tion and numerical powder averaging, including g-factor orien-
- tation effects on the dipolar coupling. Nitroxide g-tensor is
- from http://dx.doi.org/10.1063/1.1697233
- Calculation time: seconds
- Spin system properties
- Basis set
- Spinach housekeeping
- Sequence parameters
- Pulse sequence

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `state()`, `operator()`, `powder()`, `kfigure()`, `kxlabel()`.
