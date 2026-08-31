# tests/kernel/test_sparse_tensor_graph_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_sparse_tensor_graph_suite.m`
- Signature: `result=test_sparse_tensor_graph_suite()`
- Total lines: 79

## Purpose

Tests sparse, tensor-product, and simple graph utilities. Syntax: result=test_sparse_tensor_graph_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Sparse, tensor, and graph helper functions\n')`.
- Lines 20-23: State the utility target of the test; implemented by `result=new_test_result('kernel/sparse_tensor_graph_suite', 'Sparse, tensor, and graph helper functions', 'small sparse, tensor-product, graph, and combinatorial helpers…`.
- Lines 25-26: Check sparse logical matrix conversion to partial CSR indexing; implemented by `A=sparse(logical([0 1 0;1 0 1;0 0 0]))`.
- Lines 33-34: Check Kronecker-matrix multiplication without opening the product; implemented by `Q={[1 2;0 -1],[2 0;1 3],[0 1;4 -2]}`.
- Lines 42-43: Check subgraph pruning removes strict subsets while preserving maximal rows; implemented by `subgraphs=logical([1 1 0;1 0 0;0 1 1;0 1 0])`.
- Lines 48-49: Check a small permutation group table; implemented by `G=perm_group('S3')`.
- Lines 57-58: Check tuple enumeration independent of random output order; implemented by `rng(1,'twister')`.
- Lines 64-65: Check molecular connectivity on a four-point geometry; implemented by `xyz=[0 0 0;0.5 0 0;2 0 0;0 0.5 0]`.
- Lines 70-71: Check simple first-fit bin packing indices; implemented by `bins=binpack([4 2 1 5 3],5)`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/sparse_tensor_graph_suite', 'Sparse, tensor, and graph helper functions', 'small sparse, tensor-product, graph, and combinatorial helpers…`.
- Lines 26: computes `A` using `A=sparse(logical([0 1 0;1 0 1;0 0 0]))`.
- Lines 27: computes `[row_ptr,col_idx]` using `[row_ptr,col_idx]=sparse2csr(A)`.
- Lines 34: computes `Q` using `Q={[1 2;0 -1],[2 0;1 3],[0 1;4 -2]}`.
- Lines 35: computes `X` using `X=reshape(1:16,8,2)`.
- Lines 36: computes `K` using `K=kron(kron(Q{1},Q{2}),Q{3})`.
- Lines 43: computes `subgraphs` using `subgraphs=logical([1 1 0;1 0 0;0 1 1;0 1 0])`.
- Lines 44: computes `subgraphs_ref` using `subgraphs_ref=logical([1 1 0;0 1 1])`.
- Lines 49: computes `G` using `G=perm_group('S3')`.
- Lines 59: computes `tuples` using `tuples=swizzle({[1 2],[3 4 5]})`.
- Lines 60: computes `tuples_ref` using `tuples_ref=[1 3;1 4;1 5;2 3;2 4;2 5]`.
- Lines 65: computes `xyz` using `xyz=[0 0 0;0.5 0 0;2 0 0;0 0.5 0]`.
- Lines 66: computes `conn_ref` using `conn_ref=double([0 1 0 1;1 0 0 1;0 0 0 0;1 1 0 0])`.
- Lines 71: computes `bins` using `bins=binpack([4 2 1 5 3],5)`.

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
