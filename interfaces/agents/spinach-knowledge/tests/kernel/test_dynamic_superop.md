# tests/kernel/test_dynamic_superop.m

**Source:** [https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_superop.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_dynamic_superop.m)

## Purpose

Regression test for `superop()` covering sparse-XYZ spherical-tensor product superoperators. The test verifies that `superop()` produces sparse XYZ product superoperators consistent with angular-momentum algebra, checking the unit-operator shortcut, direct commutator and anticommutator identities, and Lz spherical-tensor projection eigenvalues.

## Behaviour

- Announces the test target with `fprintf('TESTING: Spherical-tensor superoperator construction\n')`.
- Initialises a test result via `new_test_result('kernel/dynamic_superop', ...)` describing the requirement that `superop()` must produce sparse XYZ product superoperators consistent with angular-momentum algebra.
- Builds a two-spin spherical-tensor Liouville-space system with `sys.magnet=0`, `sys.isotopes={'1H','13C'}`, `inter.zeeman.scalar={0,0}`, `bas.formalism='sphten-liouv'`, and `bas.approximation='none'`, using `test_spin_system(sys,inter,bas)`; the matrix dimension is taken as `size(spin_system.bas.basis,1)`.
- Unit-operator shortcut: `superop(spin_system,[0 0],'left')` is converted to sparse form and compared against `speye(matrix_dim)` with tolerances `1e-15`, asserting that an all-zero opspec must map to the unit operator over the full basis.
- For Lz on spin 1 (`[2 0]`), builds left, right, commutator (`'comm'`), and anticommutator (`'acomm'`) superoperator forms and checks:
  - `A_comm` equals `A_left - A_right` (commutator identity: a commutator superoperator must equal left multiplication minus right multiplication).
  - `A_acomm` equals `A_left + A_right` (anticommutator identity: an anticommutator superoperator must equal left multiplication plus right multiplication).
- Lz projection eigenvalues: obtains `[~,m_proj]=lin2lm(spin_system.bas.basis(:,1))` and compares `diag(A_comm)` against `m_proj`, asserting that the commutator `[Lz,T(l,m)]` must return `m*T(l,m)` on the active spin.
- Two-spin product path: for opspec `[2 2]`, builds left, right, and commutator forms and checks that `A_comm` equals `A_left - A_right` (multi-spin product commutators must obey the same sided-product identity), and that `nnz(A_comm)>0 && nnz(A_comm)<matrix_dim^2` (the two-spin product superoperator must be non-zero and sparse in the tensor basis).
- All comparisons use `test_close` with absolute and relative tolerances of `1e-15`; the sparsity check uses `test_true`.
- A local helper `local_xyz_to_sparse(xyz,matrix_dim)` converts Spinach XYZ sparse triples into MATLAB sparse storage via `A=sparse(xyz(:,1),xyz(:,2),complex(xyz(:,3)),matrix_dim,matrix_dim)`.

## Inputs and outputs

**Syntax:**

```matlab
result = test_dynamic_superop()
```

**Outputs:**

- `result` — regression test result with explanatory messages.

The function takes no inputs.

## References

- `superop()` — spherical-tensor superoperator construction under test.
- `new_test_result()`, `test_close()`, `test_true()` — test harness utilities.
- `test_spin_system()` — builds the two-spin spherical-tensor Liouville-space test system.
- `lin2lm()` — extracts the `m` projection quantum numbers from the spherical-tensor basis.
