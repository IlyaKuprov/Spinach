# examples/singlet_states/s2m_example.m

- Signature: `s2m_example()`

## Purpose

An example of the S2M sequence for a two-spin system. Calculation time: seconds

## Physical / mathematical content

- The model is a pair of 13C spins with scalar coupling 55 Hz and opposite Zeeman offsets, 0.03 and -0.03. The initial state is the two-spin singlet; the detected observable is total longitudinal magnetisation.
- The example calls the S2M sequence to convert the singlet-state preparation into a state whose longitudinal magnetisation is read out.

## Numerical / algorithmic content

- The system is represented in the sphten-liouv formalism with no basis approximation. The example constructs the NMR Hamiltonian and 13C Lx/Ly pulse operators, then calls `s2m` with the singlet initial state and parameters 55 and 6.0. It reports the overlap of the resulting state with the all-spin Lz detection state.

## Implementation structure

- An example of the S2M sequence for a two-spin system.
- Calculation time: seconds
- Spin system and interactions
- Basis set
- Spinach housekeeping
- Hamiltonian
- Pulse operators
- Start with singlet state
- Detect longitudinal magnetisation
- Call the S2M sequence
- Display the longitudinal magnetisation
