# examples/quantum_tech/transmon_duffing_ladder.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/quantum_tech/transmon_duffing_ladder.m`
- Signature: `transmon_duffing_ladder()`
- Total lines: 62

## Purpose

Duffing-model energy ladder of a weakly anharmonic transmon, showing how the transition frequencies separate as anharmo- nicity increases. Inspired by the transmon model of Koch et al., Phys. Rev. A 76, 042319 (2007). Calculation time: seconds

## Physical / mathematical content

- Quantum-technology examples. The files in this area model cavity QED, transmon qubits, NV centres, and related effective Hamiltonians. The recurring mathematics is finite-dimensional quantum dynamics with ladder operators, rotating-wave-style couplings, anharmonic oscillator terms, avoided crossings, and coherent control in coupled few-mode systems.
- The effective hardware model is a weakly anharmonic oscillator. Duffing nonlinearity breaks equal level spacing and allows qubit-like addressability within a truncated bosonic ladder.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Implementation structure

- Duffing-model energy ladder of a weakly anharmonic transmon,
- showing how the transition frequencies separate as anharmo-
- nicity increases. Inspired by the transmon model of Koch et
- al., Phys. Rev. A 76, 042319 (2007).
- Calculation time: seconds
- Magnet field
- Particle specification
- Transmon mode frequency
- Formalism and basis
- Spinach housekeeping
- Harmonic part from the declared frequency
- Anharmonicity operator

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `assume()`, `hamiltonian()`, `operator()`, `anharm()`, `trans_frq()`, `diff()`, `kfigure()`, `kxlabel()`, `kylabel()`, `ktitle()`, `klegend()`.
