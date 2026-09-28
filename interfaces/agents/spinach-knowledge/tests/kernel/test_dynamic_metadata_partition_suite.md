# tests/kernel/test_dynamic_metadata_partition_suite.m

- Signature: `result=test_dynamic_metadata_partition_suite()`

## Purpose

Regression checks for parallel-state metadata and dynamic hashing, transfer-matrix, graph-component, and disabled-helper behavior.

## Coverage

- Checks that `poolsize()` returns a nonnegative integer and `isworkernode()` is false on the MATLAB client.
- Checks stable 32-character hexadecimal `md5_hash` output, sensitivity to changed data and full versus sparse matrices, and stable first-occurrence removal of duplicate sparse rows by `unihash`.
- Checks `transfermat` recovery of a 2×2 matrix from linearly complete overdetermined samples at `1e-14` tolerance.
- Checks `scomponents` on a directed graph with components `{1,2}` and `{3}`.
- With tracking or zero-track elimination disabled, checks that `path_trace` and `zte` return the scalar placeholder `1`.

## Outputs

- `result` — regression-test result with explanatory messages.
