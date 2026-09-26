# examples/quantum_tech/circuit_qed/cross_resonance.m

- Signature: `cross_resonance()`

## Purpose

Cross-resonance gate mechanism between two fixed-frequency transmons in the laboratory frame. The control transmon is driven at the frequency of the target transmon; through the static quadrature coupling, the target then rotates about an axis in its equatorial plane by an angle conditioned on the state of the control. The difference between the two condi- tional rotations is the ZX interaction behind the native two-

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Cross-resonance gate mechanism between two fixed-frequency
- transmons in the laboratory frame. The control transmon is
- driven at the frequency of the target transmon; through the
- static quadrature coupling, the target then rotates about an
- axis in its equatorial plane by an angle conditioned on the
- state of the control. The difference between the two condi-
- tional rotations is the ZX interaction behind the native two-
- qubit gate of fixed-frequency superconducting processors; at
- these parameters the unconditional part is the larger of the
- two, and an echo sequence would be needed to remove it. A weak
- direct tone on the target sets the gate phases. The simulation
- runs in the laboratory frame, but the recorded kets are rotated
