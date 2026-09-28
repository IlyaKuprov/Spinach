# examples/fundamentals/state_tests/state_consistency_1.m

- Signature: `state_consistency_1()`

## Purpose

Test construction of two-spin singlet/triplet and four-spin product states in Spinach's Zeeman Hilbert, Zeeman Liouville, and spherical-tensor Liouville formalisms.

## Method

The script creates four `1H` spins at 14.1 T with zero Zeeman interactions and iterates over the three formalisms. For every ordered pair of distinct spins, it checks that the pair singlet plus the three triplets equals Spinach's two-spin unit state. It then checks, for every permutation of the four spins, that the 16 products formed from pair singlet/triplet states sum to the four-spin unit state. Each L1-norm residual must be at most `1e-6`; a failure raises an error, otherwise the script reports success for that formalism.
