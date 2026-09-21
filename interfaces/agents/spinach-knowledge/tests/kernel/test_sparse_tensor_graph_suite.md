# tests/kernel/test_sparse_tensor_graph_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_sparse_tensor_graph_suite.m`
- Signature: `result=test_sparse_tensor_graph_suite()`
- Total lines: 79

## Purpose

Tests sparse, tensor-product, and simple graph utilities. Syntax: result=test_sparse_tensor_graph_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `logical()`, `sparse2csr()`, `test_close()`, `kronm()`, `kronm_new()`, `test_true()`, `isequal()`, `prune_subgraphs()`, `perm_group()`, `rng()`, `swizzle()`, `sortrows()`, `double()`, `conmat()`, `binpack()`.
