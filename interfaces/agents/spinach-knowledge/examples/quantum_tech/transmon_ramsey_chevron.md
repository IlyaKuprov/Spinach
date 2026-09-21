# examples/quantum_tech/transmon_ramsey_chevron.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/transmon_ramsey_chevron.m`
- Signature: `transmon_ramsey_chevron()`
- Total lines: 86

## Purpose

Ramsey chevron of a three-level transmon in the Duffing ap- proximation. A nominal pi/2 pulse prepares a coherence, and detuning during free evolution produces Ramsey fringes. Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Implementation structure

- Ramsey chevron of a three-level transmon in the Duffing ap-
- proximation. A nominal pi/2 pulse prepares a coherence, and
- detuning during free evolution produces Ramsey fringes.
- Calculation time: seconds
- Magnet field
- Particle specification
- Transmon in the rotating frame
- Formalism and basis
- Spinach housekeeping
- Anharmonicity part from the declared interactions
- Transmon operators
- Free-evolution parameters

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `state()`, `detunings()`, `time_axis()`, `answer()`, `hdot()`, `kfigure()`, `kxlabel()`, `kylabel()`, `ktitle()`.
