# tests/kernel/test_state_constructor_suite.m

- Signature: `result=test_state_constructor_suite()`

## Purpose

Tests state-constructor helper functions.

## Physical / mathematical content

- For two spin-half nuclei, the singlet and three triplet density matrices are unit-trace projectors that sum to the identity; the singlet is idempotent and orthogonal to the zero-projection triplet.
- For two deuterons, singlet, triplet, and quintet projectors resolve the nine-dimensional Hilbert-space identity, and each population state has unit trace.
- The four-spin `S(x)S` state equals the direct product of two two-spin singlet projectors.

## Numerical / algorithmic content

- Checks one-spin unit states in Zeeman Hilbert, Zeeman Liouville, and spherical-tensor Liouville form. At finite temperature, a zero Hamiltonian gives the maximally mixed thermal state.
- Checks that `partner_state` generates four combinations for two binary partner spins and that each descriptor reproduces its corresponding `state` result.
- Checks that a zero-field triplet density matrix produced by `zftrip` has unit trace and is Hermitian.

## Outputs

- `result` — regression test result with explanatory messages.

## Implementation structure

- Announces the state-constructor test and creates a test result.
- Constructs one-, two-, three-, and four-spin test systems, a deuteron pair, and an `E3` system; compares constructor outputs with explicit states, projectors, and normalisation identities.