# kernel/utilities/add_spins.m

- Signature: `[mult,proj]=add_spins(spin_a,spin_b)`

## Purpose

Reduction of direct products of two su(2) irreps. Syntax: [mult,proj]=add_spins(spin_a,spin_b)

## Physical / mathematical content

- General mathematical and infrastructure utilities. This area contains finite differences, perturbation theory, graph algorithms, spectral densities, tensor algebra, hash/report helpers, and other reusable numerical components.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

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
