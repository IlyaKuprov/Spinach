# tests/kernel/test_sparse_tensor_graph_suite.m

- Signature: `result=test_sparse_tensor_graph_suite()`

## Purpose

Tests sparse, tensor-product, and simple graph utilities. Syntax: result=test_sparse_tensor_graph_suite()

## Physical / mathematical content

## Numerical / algorithmic content

## Outputs

- result -regression test result with explanatory messages
- The test checks small sparse-format transforms, Kronecker-matrix action,
- graph pruning, permutation group metadata, tuple enumeration, connectivity,
- and simple bin packing against explicit references.

## Implementation structure

- Tests sparse, tensor-product, and simple graph utilities. Syntax:
- result=test_sparse_tensor_graph_suite()
- result -regression test result with explanatory messages
- The test checks small sparse-format transforms, Kronecker-matrix action,
- graph pruning, permutation group metadata, tuple enumeration, connectivity,
- and simple bin packing against explicit references.
- Announce the test target
- State the utility target of the test
- Check sparse logical matrix conversion to partial CSR indexing
- Check Kronecker-matrix multiplication without opening the product
- Check subgraph pruning removes strict subsets while preserving maximal rows
- Check a small permutation group table
