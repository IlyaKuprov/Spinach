# tests/kernel/test_sparse_tensor_utility_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_sparse_tensor_utility_suite.m`
- Signature: `result=test_sparse_tensor_utility_suite()`
- Total lines: 69

## Purpose

Tests sparse, tensor, and numerical utility helpers. Syntax: result=test_sparse_tensor_utility_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Sparse, tensor, and numerical utility functions\n')`.
- Lines 20-23: State the utility target of the test; implemented by `result=new_test_result('kernel/sparse_tensor_utility_suite', 'Sparse, tensor, and numerical utility functions', 'sparse/tensor utilities must match direct dense algebra…`.
- Lines 25-26: Kronecker-product application must match explicit kron products; implemented by `Q={sparse([1 2;3 4]),sparse([0 5;6 7])}`.
- Lines 34-35: Blicharski invariants ignore isotropic trace and follow their definitions; implemented by `A=[1 2 3;4 5 6;7 8 9]`.
- Lines 52-53: Lorentzian spectral density follows its closed definition; implemented by `L=2; Drot=1.5e6; omega=2.0e5; tau=1/(L*(L+1)*Drot)`.
- Lines 57-58: Frobenius SVD truncation keeps the smallest rank whose dropped tail is below tolerance; implemented by `s=[5 1 0.01]`.
- Lines 62-63: svd_shrink factorisation must reconstruct retained singular content; implemented by `spin_system.sys.output='hush'`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/sparse_tensor_utility_suite', 'Sparse, tensor, and numerical utility functions', 'sparse/tensor utilities must match direct dense algebra…`.
- Lines 26: computes `Q` using `Q={sparse([1 2;3 4]),sparse([0 5;6 7])}`.
- Lines 27: computes `x` using `x=(1:4)'`.
- Lines 28: computes `K` using `K=kron(Q{1},Q{2})`.
- Lines 35: computes `A` using `A=[1 2 3;4 5 6;7 8 9]`.
- Lines 36: computes `[Lsq,Dsq]` using `[Lsq,Dsq]=blinv(A)`.
- Lines 37: computes `Lsq_ref` using `Lsq_ref=(A(1,2)-A(2,1))^2+(A(1,3)-A(3,1))^2+(A(2,3)-A(3,2))^2`.
- Lines 38-39: computes `Dsq_ref` using `Dsq_ref=A(1,1)^2+A(2,2)^2+A(3,3)^2-A(1,1)*A(2,2)-A(1,1)*A(3,3)-A(2,2)*A(3,3)+ (3/4)*((A(1,2)+A(2,1))^2+(A(1,3)+A(3,1))^2+(A(2,3)+A(3,2))^2)`.
- Lines 44: computes `B` using `B=[2 -1 0; 3 4 1; 5 6 -2]`.
- Lines 45: computes `[X1,X2]` using `[X1,X2]=blprod(A,B)`.
- Lines 46: computes `[Lam,Dam]` using `[Lam,Dam]=blinv(A-B); [Lap,Dap]=blinv(A+B)`.
- Lines 53: computes `L` using `L=2; Drot=1.5e6; omega=2.0e5; tau=1/(L*(L+1)*Drot)`.
- Lines 58: computes `s` using `s=[5 1 0.01]`.
- Lines 63: computes `spin_system.sys.output` using `spin_system.sys.output='hush'`.
- Lines 64: computes `rho` using `rho=diag([4 1 1e-6])`.
- Lines 65: computes `[vec,cov]` using `[vec,cov]=svd_shrink(spin_system,rho,1e-4)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks Kronecker-product application, Blicharski invariants,
- spectral densities, SVD truncation helpers, sparse density, and related
- numerical utilities against direct definitions.

## Implementation structure

- Tests sparse, tensor, and numerical utility helpers. Syntax:
- result=test_sparse_tensor_utility_suite()
- result -regression test result with explanatory messages
- The test checks Kronecker-product application, Blicharski invariants,
- spectral densities, SVD truncation helpers, sparse density, and related
- numerical utilities against direct definitions.
- Announce the test target
- State the utility target of the test
- Kronecker-product application must match explicit kron products
- Blicharski invariants ignore isotropic trace and follow their definitions
- Lorentzian spectral density follows its closed definition
- Frobenius SVD truncation keeps the smallest rank whose dropped tail is below tolerance

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `test_close()`, `kronm()`, `kronm_new()`, `blinv()`, `blprod()`, `spden()`, `frob_chop()`, `svd_shrink()`.
