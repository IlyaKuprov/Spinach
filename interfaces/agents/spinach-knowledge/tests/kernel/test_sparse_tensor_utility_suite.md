# tests/kernel/test_sparse_tensor_utility_suite.m

## Purpose

Regression test suite for the sparse, tensor, and numerical utility helpers in Spinach. The suite verifies that Kronecker-product application, Blicharski invariants, spectral densities, SVD truncation helpers, sparse density, and related numerical utilities agree with their direct dense-algebra definitions on small reference cases.

Source: [tests/kernel/test_sparse_tensor_utility_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_sparse_tensor_utility_suite.m)

## Behaviour

The function announces the test target with `fprintf`, initialises a test result object via `new_test_result` under the identifier `kernel/sparse_tensor_utility_suite`, and then runs a sequence of `test_close` comparisons:

- **Kronecker-product application**: builds two sparse 2×2 matrices `Q{1}` and `Q{2}` and a 4-element vector `x=(1:4)'`, forms the explicit product `K=kron(Q{1},Q{2})`, and checks that `kronm(Q,x)` and `kronm_new(Q,x)` both match `K*x` to absolute and relative tolerances of `1e-14`.
- **Blicharski invariants**: for `A=[1 2 3;4 5 6;7 8 9]`, checks that `blinv` returns the first-rank invariant `Lsq` equal to the sum of squared antisymmetric off-diagonal differences, and the second-rank invariant `Dsq` equal to the traceless symmetric tensor amplitude (diagonal combination plus `3/4` times the sum of squared symmetric off-diagonal sums). For a second matrix `B=[2 -1 0; 3 4 1; 5 6 -2]`, checks that `blprod(A,B)` returns first- and second-rank tensor products `X1` and `X2` matching the polarisation identities `(Lap-Lam)/4` and `(Dap-Dam)/4`, where `blinv` is evaluated at `A-B` and `A+B`. All comparisons use tolerances of `1e-14`.
- **Spectral density**: with `L=2`, `Drot=1.5e6`, `omega=2.0e5`, and `tau=1/(L*(L+1)*Drot)`, checks that `spden(L,Drot,omega)` equals `(tau/(2*L+1))/(1+(tau*omega)^2)` with absolute tolerance `1e-20` and relative tolerance `1e-14`.
- **Frobenius SVD truncation**: for the singular value vector `s=[5 1 0.01]` and tolerance `0.02`, checks that `frob_chop(s,0.02)` returns `2` exactly (zero tolerances), i.e. only the `0.01` singular value may be dropped.
- **SVD shrink factorisation**: with `spin_system.sys.output='hush'` and `rho=diag([4 1 1e-6])`, calls `svd_shrink(spin_system,rho,1e-4)` and checks that `vec*cov'` reconstructs `diag([4 1 0])` to tolerances of `1e-12`, i.e. singular values above the threshold are retained as vector/covector factors.

## Inputs and outputs

**Syntax**

```matlab
result = test_sparse_tensor_utility_suite()
```

The function takes no inputs.

**Outputs**

- `result` — regression test result object with explanatory messages, accumulated through the `test_close` comparisons described above.

## References

- [Spinach GitHub repository — tests/kernel/test_sparse_tensor_utility_suite.m](https://github.com/IlyaKuprov/Spinach/blob/main/tests/kernel/test_sparse_tensor_utility_suite.m)
- [Spinach project website](https://spindynamics.org/)
