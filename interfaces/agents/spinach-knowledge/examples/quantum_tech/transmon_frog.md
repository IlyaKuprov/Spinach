# examples/quantum_tech/transmon_frog.m

- Signature: `transmon_frog()`

## Purpose

Basic implementation of a Frequency Robust Gate (FROG) for a single transmon, based on: Ensemble GRAPE optimisation with a distribution of control powers and excess amplitude penalty. Calculation time: minutes

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.
- Orientation or trajectory averaging is performed numerically, so grid design, weights, and integration error control matter directly to accuracy and runtime.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Basic implementation of a Frequency Robust Gate
- (FROG) for a single transmon, based on:
- Ensemble GRAPE optimisation with a distribution
- of control powers and excess amplitude penalty.
- Calculation time: minutes
- Magnet field
- Particle specification
- Rotating frame transmon parameters
- Formalism and basis
- Spinach housekeeping
- Drift Hamiltonian from the declared interactions
- Build the intial control pulses (FROG)
