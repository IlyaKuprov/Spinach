# examples/quantum_tech/circuit_qed/transmon_drag.m

- Signature: `transmon_drag()`

## Purpose

DRAG correction of a resonant Gaussian pulse on a three-level Duffing transmon in the laboratory frame. The derivative of the envelope, scaled by the DRAG detuning parameter, is fed into the quadrature channel of the IQ mixer; this suppresses the leakage into the second excited state and improves the fidelity of the 90 degree rotation on the qubit subspace. A further improvement comes from numerical optimisation of t

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- DRAG correction of a resonant Gaussian pulse on a three-level
- Duffing transmon in the laboratory frame. The derivative of the
- envelope, scaled by the DRAG detuning parameter, is fed into the
- quadrature channel of the IQ mixer; this suppresses the leakage
- into the second excited state and improves the fidelity of the
- 90 degree rotation on the qubit subspace. A further improvement
- comes from numerical optimisation of the mixer parameters. The
- rows of pulse_params are the plain Gaussian, the analytically
- corrected DRAG, and the numerically optimised DRAG pulse; the
- columns are the envelope amplitude, the DRAG detuning, the lo-
- cal oscillator frequency, the mixer phase, and the DRAG quad-
- rature switch. The propagator is taken into the frame of the
