# tests/kernel/test_sparse_tensor_graph_suite.m

- Signature: `result=test_sparse_tensor_graph_suite()`

## Purpose

Checks sparse-format conversion, Kronecker-product operations, subgraph pruning, permutation-group metadata, tuple enumeration, connectivity, and bin packing.

## Physical / mathematical content

Not applicable; this is a utility regression suite.

## Numerical / algorithmic content

The test converts a 3-by-3 sparse logical matrix to one-based CSR row pointers and indices, and checks `kronm` and `kronm_new` against explicit Kronecker-product action on an 8-by-2 matrix. It checks that `prune_subgraphs` keeps maximal rows, verifies a small permutation-group table, and compares enumeration of `[1 2]` with `[3 4 5]` to the six expected tuples. It also checks the connectivity matrix for four points at Euclidean cutoff `0.75`, and verifies greedy bin packing of sizes `[4 2 1 5 3]` at capacity `5` into four bins.

## Outputs

`result` is the regression-test record with explanatory messages.

## Implementation structure

Runs each utility on a small explicit fixture and compares its output with the expected indices, matrices, tuples, or bins.
