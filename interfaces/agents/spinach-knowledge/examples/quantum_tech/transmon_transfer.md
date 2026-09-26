# examples/quantum_tech/transmon_transfer.m

- Signature: `transmon_transfer()`

## Purpose

Basic two-transmon system with Duffing model interacti- ons and a flip-flop coupling; coherence transfer from transmon 1 to transmon 2. GRAPE optimisation with a distribution control powers and transmon offsets with a penalty on excess power. Calculation time: minutes.

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- Basic two-transmon system with Duffing model interacti-
- ons and a flip-flop coupling; coherence transfer from
- transmon 1 to transmon 2. GRAPE optimisation with a
- distribution control powers and transmon offsets with
- a penalty on excess power.
- Calculation time: minutes.
- Magnet field
- Particle specification
- Rotating frame transmon parameters
- Formalism and basis
- Spinach housekeeping
- Drift Hamiltonian from the declared interactions
