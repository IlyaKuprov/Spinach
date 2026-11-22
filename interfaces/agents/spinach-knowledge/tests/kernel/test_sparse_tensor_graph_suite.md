# tests/kernel/test_sparse_tensor_graph_suite.m

**Source:** https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_sparse_tensor_graph_suite.m

## Purpose

Regression test suite for the small sparse-format, tensor-product, graph, and combinatorial helper functions in the kernel module. The suite verifies deterministic reference behaviour of these helpers against explicit references.

## Behaviour

The function announces the test target with `fprintf`, initialises a result object via `new_test_result` for the target `kernel/sparse_tensor_graph_suite`, and then runs a sequence of checks:

- **`sparse2csr`** — converts a sparse logical 3-by-3 matrix to partial CSR indexing and checks that the row pointer equals `[1;2;4;4]` (one-based row starts plus a final sentinel) and the column index equals `[2;1;3]` (row-major non-zero listing), each with tolerances `1e-15`.
- **`kronm` and `kronm_new`** — for the cell array `Q={[1 2;0 -1],[2 0;1 3],[0 1;4 -2]}` and the 8-by-2 matrix `X=reshape(1:16,8,2)`, both functions are checked against multiplication by the explicit triple Kronecker product `kron(kron(Q{1},Q{2}),Q{3})`, with tolerances `1e-14`.
- **`prune_subgraphs`** — from the logical matrix `[1 1 0;1 0 0;0 1 1;0 1 0]`, strict subsets are removed while maximal rows are preserved, giving the reference `[1 1 0;0 1 1]`.
- **`perm_group('S3')`** — checks class sizes `[1 2 3]` (identity, two 3-cycles, three transpositions), the irreducible character table `[1 1 1;1 1 -1;2 -1 0]`, and that the group order equals the sum of its conjugacy class sizes.
- **`swizzle`** — with `rng(1,'twister')` set beforehand, enumerates the Cartesian product of `{[1 2],[3 4 5]}`; because output order may be random, the result is sorted with `sortrows` and compared to `[1 3;1 4;1 5;2 3;2 4;2 5]` with tolerances `1e-15`.
- **`conmat`** — on the four-point geometry `xyz=[0 0 0;0.5 0 0;2 0 0;0 0.5 0]` with cutoff `0.75`, checks the connectivity matrix against the reference `[0 1 0 1;1 0 0 1;0 0 0 0;1 1 0 0]`, i.e. points whose Euclidean separation is below the cutoff are connected.
- **`binpack`** — for the input `[4 2 1 5 3]` with bin capacity `5`, checks that four bins are returned with contents `1`, `[2;3]`, `4`, and `5` respectively, confirming greedy first-fit filling from the remaining list with original indices returned.

Each check appends an explanatory message to the test result via `test_close` or `test_true`.

## Inputs and outputs

**Syntax:**

```matlab
result = test_sparse_tensor_graph_suite()
```

**Outputs:**

- `result` — regression test result with explanatory messages.

The function takes no inputs.

## References

- Source file: https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_sparse_tensor_graph_suite.m
