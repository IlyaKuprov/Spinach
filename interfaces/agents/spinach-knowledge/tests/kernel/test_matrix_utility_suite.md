# tests/kernel/test_matrix_utility_suite.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_matrix_utility_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_matrix_utility_suite.m)

## Purpose

Regression test suite for the small matrix utility functions that Spinach uses throughout the kernel for algebra, sparsity, block assembly, and indexing operations. The suite verifies that these low-level matrix helpers preserve their exact algebraic definitions.

## Behaviour

The function announces the test target with `fprintf('TESTING: Matrix utility functions\n')` and initialises a test result object via `new_test_result('kernel/matrix_utility_suite', 'Matrix utility functions', 'matrix helpers must preserve their exact algebraic definitions.')`. It then runs a sequence of checks:

- **Anticommutator:** `acomm(A,B)` is compared against `A*B+B*A` for the test matrices `A=[1 2;3 4]` and `B=[0 1;-1 2]`.
- **Cheap norm (CPU):** `cheap_norm(A)` is compared against `norm(A,1)` for a CPU matrix.
- **Polyadic norm estimation:** `cheap_norm` is applied to `polyadic({{A}})` and compared against `norm(A,1)`; the block estimator `cheap_norm(P,2,5)` is also checked. Tall (`R=[1 0;2 3;4 5]`) and wide (`R=[1 2 3;4 5 6]`) rectangular polyadics are checked to use row-length sign vectors and column-length probes.
- **Higham–Tisseur zero-sign convention:** for `Z=[-2 0;-1 2]`, zero phase entries are treated as one, as required by Algorithm 2.4.
- **Complex adjoint and phase handling:** for `Z=[-2 0;-1+1i 2]`, complex polyadic estimation uses conjugate transposes and unit phases.
- **Block estimator lower-bound contract:** with `rng(1)` (state saved and restored via `onCleanup`), for `Z=[0 -3 2 1;4 0 -1 2;-2 5 0 -4;1 -1 3 0]`, `cheap_norm(P,3,5)` must return a positive lower bound not exceeding the exact one-norm (checked as `est<=norm(Z,1)+1e-12` and `est>0`).
- **Polyadic realness dispatch:** `isreal(polyadic({{Z}}))` is true for real cores and false for a complex core `[1 1i;0 2]`, recognised without opening the Kronecker products.
- **Identity and trace predicates:** `iseye(speye(3))` is true; `iseye([1 0;0 2])` is false; `istraceless([1 0;0 -1])` is true; `istraceless(eye(2))` is false; `krondelta(3,3)` is true; `krondelta(2,3)` is false.
- **Sparse block-diagonal assembly:** `sp_block_diag([1 2;3 4],[5;6])` is compared against `sparse([1 2 0;3 4 0;0 0 5;0 0 6])`.
- **Row and column replication:** `repcols(M,2,3)` for `M=[1 2 3;4 5 6]` is compared against `[1 2 2 2 3;4 5 5 5 6]`; `reprows(M,1,2)` is compared against `[1 2 3;1 2 3;4 5 6]`.
- **Sparse cleanup:** with a spin system configured with `sys.output='hush'`, `sys.disable={}`, `sys.enable={}`, `tols.dense_matrix=0.9`, `tols.small_matrix=10`, `clean_up(spin_system,C,1e-10)` applied to `C=sparse([1e-12 1;0 2])` is compared against `sparse([0 1;0 2])`, verifying that elements below the requested non-zero tolerance are removed.

All closeness checks use tolerances of `1e-15` for both the absolute and relative tolerances.

## Inputs and outputs

- **Outputs:**
  - `result` — regression test result with explanatory messages.
- **Inputs:** none.

## References

- [Spinach MATLAB source on GitHub](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_matrix_utility_suite.m)
