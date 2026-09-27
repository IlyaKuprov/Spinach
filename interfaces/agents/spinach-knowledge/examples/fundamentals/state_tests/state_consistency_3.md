# examples/fundamentals/state_tests/state_consistency_3.m

- Signature: `state_consistency_3()`

## Purpose

Check the deuterium-pair singlet, triplet, and quintet projectors returned by `deut_pair`, and verify their completeness.

## Method

For a pair of `2H` spins with zero field and zero Zeeman interactions, the script compares the Hilbert-space singlet, three triplet, and five quintet projectors against explicit component-state projectors. The component definitions follow Eq. 1 in [the cited paper](https://doi.org/10.1016/S0009-2614(98)00784-2). The source then rebuilds the system with the `zeeman-hilb` formalism and checks that the singlet, all three triplets, and all five quintets sum to the two-spin unit state. Residuals above `1e-6` fail the test.
