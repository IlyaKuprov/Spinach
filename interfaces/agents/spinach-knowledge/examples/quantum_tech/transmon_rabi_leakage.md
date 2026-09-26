# examples/quantum_tech/transmon_rabi_leakage.m

- Signature: `transmon_rabi_leakage()`

## Purpose

Rabi dynamics of a driven four-level transmon in the Duffing approximation, including leakage into the second and third excited states. The resonant drive is a part of the rotating frame Hamiltonian, and all four level populations come from a single trajectory. Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

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
