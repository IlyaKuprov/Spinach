# examples/quantum_tech/transmon_rabi_leakage.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/transmon_rabi_leakage.m`
- Signature: `transmon_rabi_leakage()`
- Total lines: 61

## Purpose

Rabi dynamics of a driven four-level transmon in the Duffing approximation, including leakage into the second and third excited states. The resonant drive is a part of the rotating frame Hamiltonian, and all four level populations come from a single trajectory. Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 13-14: Magnet field; implemented by `sys.magnet=0`.
- Lines 16-17: Particle specification; implemented by `sys.isotopes={'T4'}`.
- Lines 19-20: Resonantly driven transmon in the rotating frame; implemented by `inter.modes.frqs={0}`.
- Lines 23-24: Formalism and basis; implemented by `bas.formalism='zeeman-hilb'`.
- Lines 27-28: Spinach housekeeping; implemented by `spin_system=create(sys,inter)`.
- Lines 31-32: Drift Hamiltonian from the declared interactions; implemented by `H=hamiltonian(assume(spin_system,'cavity'))`.
- Lines 34-35: Resonant drive on the transmon quadrature; implemented by `Cr=operator(spin_system,'C',1)`.
- Lines 39-40: Initial condition; implemented by `rho=state(spin_system,'BL1',1)`.
- Lines 42-43: Trajectory of the driven transmon; implemented by `traj=evolution(spin_system,H,[],rho,1e-9,400,'trajectory')`.
- Lines 45-46: Level populations; implemented by `pops=zeros(4,401)`.
- Lines 52-53: Plot the leakage dynamics; implemented by `time_axis=linspace(0,400,401)`.

### Control flow inferred from the code

- Line 47: `for` loop over `n=1:4`.

### Key state/data transformations

- Lines 14: computes `sys.magnet` using `sys.magnet=0`.
- Lines 17: computes `sys.isotopes` using `sys.isotopes={'T4'}`.
- Lines 20: computes `inter.modes.frqs` using `inter.modes.frqs={0}`.
- Lines 21: computes `inter.modes.anharms` using `inter.modes.anharms={-250e6}`.
- Lines 24: computes `bas.formalism` using `bas.formalism='zeeman-hilb'`.
- Lines 25: computes `bas.approximation` using `bas.approximation='none'`.
- Lines 28: computes `spin_system` using `spin_system=create(sys,inter)`.
- Lines 32: computes `H` using `H=hamiltonian(assume(spin_system,'cavity'))`.
- Lines 35: computes `Cr` using `Cr=operator(spin_system,'C',1)`.
- Lines 36: computes `An` using `An=operator(spin_system,'A',1)`.
- Lines 40: computes `rho` using `rho=state(spin_system,'BL1',1)`.
- Lines 43: computes `traj` using `traj=evolution(spin_system,H,[],rho,1e-9,400,'trajectory')`.
- Lines 46: computes `pops` using `pops=zeros(4,401)`.
- Lines 48: computes `coil` using `coil=state(spin_system,['BL' int2str(n)],1)`.
- Lines 49: computes `pops(n,:)` using `pops(n,:)=real(cellfun(@(rho)full(hdot(coil,rho)),traj))`.
- Lines 53: computes `time_axis` using `time_axis=linspace(0,400,401)`.

## Implementation structure

- Rabi dynamics of a driven four-level transmon in the Duffing
- approximation, including leakage into the second and third
- excited states. The resonant drive is a part of the rotating
- frame Hamiltonian, and all four level populations come from
- a single trajectory.
- Calculation time: seconds
- Magnet field
- Particle specification
- Resonantly driven transmon in the rotating frame
- Formalism and basis
- Spinach housekeeping
- Drift Hamiltonian from the declared interactions

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `state()`, `evolution()`, `int2str()`, `pops()`, `cellfun()`, `hdot()`, `kfigure()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
