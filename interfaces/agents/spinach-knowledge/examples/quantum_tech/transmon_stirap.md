# examples/quantum_tech/transmon_stirap.m

- Signature: `transmon_stirap()`

## Purpose

Basic single transmon system with Duffing model in- teractions, parameters and model from: Ensemble GRAPE optimisation with a distribution of control powers and a penalty on excess power. Calculation time: minutes

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Basic single transmon system with Duffing model in-
- teractions, parameters and model from:
- Ensemble GRAPE optimisation with a distribution of
- control powers and a penalty on excess power.
- Calculation time: minutes
- Magnet field
- Particle specification
- Rotating frame ladder detunings
- Formalism and basis
- Spinach housekeeping
- Drift Hamiltonian from the declared interactions
- Pulse power ensemble
