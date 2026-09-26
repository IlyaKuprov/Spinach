# examples/quantum_tech/circuit_qed/transmon_two_photon.m

- Signature: `transmon_two_photon()`

## Purpose

Two-photon transition in a four-level Duffing transmon. The drive carrier is placed halfway between the 0-1 and 1-2 tran- sition frequencies, where neither single-photon transition is resonant, and GRAPE finds a pulse that moves the population from the ground state into the second excited state through a virtual intermediate state. Model and parameters from Example 2 of the GRAPE_SCQ package. Calculation time: minute

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The numerical method is quasi-Newton optimisation: curvature information is approximated from successive step and gradient differences instead of forming exact second derivatives every iteration.
- The numerical method is limited-memory quasi-Newton optimisation, which keeps only a short curvature history and is therefore suitable for waveform vectors too large for dense Hessians.
- The optimisation logic is Newton or Newton-like: search directions use first- and second-order local curvature information, usually with regularisation or line-search safeguards.
- The control theory content is GRAPE: fidelity derivatives are propagated through a piecewise-constant pulse sequence so that waveform samples can be improved by gradient-based optimisation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.
- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Two-photon transition in a four-level Duffing transmon. The
- drive carrier is placed halfway between the 0-1 and 1-2 tran-
- sition frequencies, where neither single-photon transition is
- resonant, and GRAPE finds a pulse that moves the population
- from the ground state into the second excited state through a
- virtual intermediate state. Model and parameters from Example
- 2 of the GRAPE_SCQ package.
- Calculation time: minutes
- Magnet field
- Particle specification
- Transmon detuning from the carrier and anharmonicity
- Formalism and basis
