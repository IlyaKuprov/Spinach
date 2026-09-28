# tests/kernel/test_giant_ham_descr.m

- Signature: `result=test_giant_ham_descr()`

## Purpose

Checks that the giant-spin Hamiltonian descriptor route agrees with direct spherical-tensor assembly for high-rank terms.

## Physical / mathematical content

The test builds a compact giant-spin system with tensor coefficients and Euler angles, then compares the descriptor-generated Hamiltonian with a direct spherical-tensor reference after orientation. It exercises complete and secular giant-spin terms and a separate assumption-specific case.

## Numerical / algorithmic content

The local test helper constructs the spin system, applies the requested Hamiltonian assumption, assembles the production Hamiltonian and orientation contribution, and compares it with the direct tensor sum. The spherical-tensor projection loop contracts coefficient components and adds each operator contribution.

## Outputs

`result` is the regression-test result with explanatory messages.

## Implementation structure

The top-level test calls a local case helper for complete and secular giant-spin terms and for the `deer-zz` secular case. The helper evaluates the descriptor route and independently summed spherical-tensor contributions, including the Hermitisation convention used by orientation assembly.
