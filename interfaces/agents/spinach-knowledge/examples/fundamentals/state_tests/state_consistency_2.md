# examples/fundamentals/state_tests/state_consistency_2.m

- Signature: `state_consistency_2()`

## Purpose

Verify that a constructed two-spin state is represented consistently in Zeeman Hilbert space, Zeeman Liouville space, and spherical-tensor Liouville space after projection back to the Zeeman basis.

## Method

The example uses non-interacting `14N` and `235U` spins at 14.1 T. In each formalism it builds the same operator state: Lz on spin 1 and Lx on spin 2, plus L+ on spin 1. The Zeeman Liouville representation is reshaped to a 24-by-24 matrix; the spherical-tensor representation is first transformed with `sphten2zeeman` and then reshaped. The three matrices must agree pairwise within an L1-norm tolerance of `1e-6`, or the test fails.
