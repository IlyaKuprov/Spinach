# examples/esr_sol_pulsed/hard_3_pulse_deer_exchange.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/esr_sol_pulsed/hard_3_pulse_deer_exchange.m`
- Signature: `hard_3_pulse_deer_exchange()`
- Total lines: 88

## Purpose

Three-pulse DEER on a Cu(II)-Cu(II) system in a linked porphyrin complex with a strong exchange coupling between the electrons. A distribution in the exchange coupling is summed over. The calculation is done by brute-force time propagation and numerical powder averaging in Liouville space. Calculation time: minutes

## Physical / mathematical content

- Pulsed ESR / EPR solid-state examples. These scripts revolve around electron spin echo sequences, DEER, RIDME, ENDOR, ESEEM, and HYSCORE. They combine anisotropic Zeeman and hyperfine Hamiltonians with selective pulses, echo formation, and orientation averaging.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Generate the distribution; implemented by `j_values=linspace(6e6,20e6,20)`.
- Lines 21-22: Run the averaging; implemented by `for n=1:numel(weights)`.
- Lines 24-25: Hush up; implemented by `sys.output='hush'`.
- Lines 27-28: Magnet field; implemented by `sys.magnet=1.2132`.
- Lines 30-31: Isotopes; implemented by `sys.isotopes={'E','E'}`.
- Lines 33-34: Zeeman interactions; implemented by `inter.zeeman.eigs=cell(1,2)`.
- Lines 41-42: Exchange coupling; implemented by `inter.coupling.scalar{1,2}=j_values(n)`.
- Lines 45-46: Coordinates; implemented by `inter.coordinates=cell(2,1)`.
- Lines 50-51: Basis set; implemented by `bas.formalism='sphten-liouv'`.
- Lines 54-55: Disable trajectory level SSR algorithms; implemented by `sys.disable={'hygiene','trajlevel'}`.
- Lines 57-58: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 61-62: Sequence parameters; implemented by `parameters.rho0=state(spin_system,'Lz','E')`.
- Lines 72-73: Pulse sequence; implemented by `deer=powder(spin_system,@deer_3p_hard_deer,parameters,'deer-zz')`.
- Lines 75-76: Add to the total; implemented by `answer=answer+weights(n)*deer.deer_trace`.
- Lines 80-81: Build time axis; implemented by `time_axis=linspace(0,parameters.stepsize*parameters.nsteps,parameters.nsteps+1)`.
- Lines 83-84: Plotting; implemented by `kfigure(); plot(1e6*time_axis,imag(answer))`.

### Control flow inferred from the code

- Line 22: `for` loop over `n=1:numel(weights)`.

### Key state/data transformations

- Lines 17: computes `j_values` using `j_values=linspace(6e6,20e6,20)`.
- Lines 18: computes `weights` using `weights=gaussfun(j_values-13.1e6,4.2622e6)`.
- Lines 25: computes `sys.output` using `sys.output='hush'`.
- Lines 28: computes `sys.magnet` using `sys.magnet=1.2132`.
- Lines 31: computes `sys.isotopes` using `sys.isotopes={'E','E'}`.
- Lines 34: computes `inter.zeeman.eigs` using `inter.zeeman.eigs=cell(1,2)`.
- Lines 35: computes `inter.zeeman.euler` using `inter.zeeman.euler=cell(1,2)`.
- Lines 36: computes `inter.zeeman.eigs{1}` using `inter.zeeman.eigs{1}=[2.050 2.050 2.195]`.
- Lines 37: computes `inter.zeeman.euler{1}` using `inter.zeeman.euler{1}=[0 0 0]`.
- Lines 38: computes `inter.zeeman.eigs{2}` using `inter.zeeman.eigs{2}=[2.050 2.050 2.195]`.
- Lines 39: computes `inter.zeeman.euler{2}` using `inter.zeeman.euler{2}=[0 0 0]`.
- Lines 42: computes `inter.coupling.scalar{1,2}` using `inter.coupling.scalar{1,2}=j_values(n)`.
- Lines 43: computes `inter.coupling.scalar{2,2}` using `inter.coupling.scalar{2,2}=[]`.
- Lines 46: computes `inter.coordinates` using `inter.coordinates=cell(2,1)`.
- Lines 47: computes `inter.coordinates{1}` using `inter.coordinates{1}=[0.00 0.00 0.00]`.
- Lines 48: computes `inter.coordinates{2}` using `inter.coordinates{2}=[24.50 0.00 0.00]`.
- Lines 51: computes `bas.formalism` using `bas.formalism='sphten-liouv'`.
- Lines 52: computes `bas.approximation` using `bas.approximation='none'`.

## Implementation structure

- Three-pulse DEER on a Cu(II)-Cu(II) system in a linked porphyrin
- complex with a strong exchange coupling between the electrons. A
- distribution in the exchange coupling is summed over.
- The calculation is done by brute-force time propagation and
- numerical powder averaging in Liouville space.
- Calculation time: minutes
- Generate the distribution
- Run the averaging
- Hush up
- Magnet field
- Isotopes
- Zeeman interactions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `gaussfun()`, `j_values()`, `create()`, `basis()`, `state()`, `operator()`, `powder()`, `weights()`, `kfigure()`, `kxlabel()`.
