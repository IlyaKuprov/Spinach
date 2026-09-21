# kernel/utilities/add_spins.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/utilities/add_spins.m`
- Signature: `[mult,proj]=add_spins(spin_a,spin_b)`
- Total lines: 106

## Purpose

Reduction of direct products of two su(2) irreps. Syntax: [mult,proj]=add_spins(spin_a,spin_b)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_a -quantum number of the first spin,
- an integer or a half-integer
- spin_b -quantum number of the second spin,
- an integer or a half-integer

## Outputs

- mult -multiplicities corresponding to
- the values of the total spin that
- are present
- proj -projectors that reduce the direct
- product representation, a cell ar-
- ray of matrices

## Implementation structure

- Reduction of direct products of two su(2) irreps. Syntax:
- [mult,proj]=add_spins(spin_a,spin_b)
- spin_a -quantum number of the first spin,
- an integer or a half-integer
- spin_b -quantum number of the second spin,
- mult -multiplicities corresponding to
- the values of the total spin that
- are present
- proj -projectors that reduce the direct
- product representation, a cell ar-
- ray of matrices
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `pauli()`, `uint32()`, `mult()`, `any()`, `Sx_irr()`, `isscalar()`.
