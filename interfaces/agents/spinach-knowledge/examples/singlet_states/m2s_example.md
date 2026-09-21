# examples/singlet_states/m2s_example.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/singlet_states/m2s_example.m`
- Signature: `m2s_example()`
- Total lines: 45

## Purpose

An example of the M2S sequence for a two-spin system. Calculation time: seconds

## Physical / mathematical content

- Long-lived singlet-state examples. The central concept is symmetry-protected or nearly symmetry-protected two-spin order that relaxes much more slowly than ordinary Zeeman magnetisation. Files here often analyse singlet-triplet subspaces, state conversion sequences, and relaxation leakage channels.
- The relevant state manifold is the singlet/triplet decomposition, where permutation symmetry controls selection rules, relaxation susceptibility, and convertibility to ordinary magnetisation.

## Numerical / algorithmic content

- The file is built around the standard Spinach workflow: create the spin system, choose a basis or context, assemble operators/superoperators, then propagate or analyse the resulting dynamics.

## Implementation structure

- An example of the M2S sequence for a two-spin system.
- Calculation time: seconds
- Spin system and interactions
- Basis set
- Spinach housekeeping
- Hamiltonian
- Pulse operators
- Start with longitudinal magnetisation
- Detect singlet state
- Call the M2S sequence
- Display the singlet population

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `create()`, `basis()`, `hamiltonian()`, `assume()`, `operator()`, `state()`, `singlet()`, `m2s()`, `num2str()`.
