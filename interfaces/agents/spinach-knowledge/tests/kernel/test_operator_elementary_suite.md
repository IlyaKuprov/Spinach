# tests/kernel/test_operator_elementary_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_operator_elementary_suite.m`
- Signature: `result=test_operator_elementary_suite()`
- Total lines: 103

## Purpose

Tests elementary operator generators in kernel/operators. Syntax: result=test_operator_elementary_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 16-17: Announce the test target; implemented by `fprintf('TESTING: Elementary operator generators\n')`.
- Lines 19-22: State the operator target of the test; implemented by `result=new_test_result('kernel/operator_elementary_suite', 'Elementary operator generators', 'low-level operator constructors must satisfy their defining algebraic ident…`.
- Lines 24-25: Check spin-one angular momentum commutation and ladder definitions; implemented by `S=pauli(3)`.
- Lines 33-34: Check finite-truncation Weyl algebra away from the unavoidable edge state; implemented by `W=weyl(4)`.
- Lines 42-43: Check bosonic monomial serpentine indexing; implemented by `B=boson_mono(3)`.
- Lines 51-52: Check Gram-Schmidt orthogonality without imposing normalisation; implemented by `B=boson_ortho(3)`.
- Lines 62-63: Check single-transition basis indexing from the documented 4x4 map; implemented by `A=sin_tran(4)`.
- Lines 71-72: Check central-transition operators embedded in a spin-3/2 manifold; implemented by `ct_z=spalloc(4,4,2); ct_z(2,2)=0.5; ct_z(3,3)=-0.5`.
- Lines 82-83: Check irreducible spherical tensor projection quantum numbers; implemented by `T=irr_sph_ten(3,2)`.
- Lines 91-93: Check the simplest Stevens operator and Hermiticity of observable components; implemented by `result=test_close(result,'stevens rank zero',stevens(3,0,0),speye(3),1e-15,1e-15, 'the rank-zero Stevens operator is the unit operator')`.

### Control flow inferred from the code

- Line 54: `for` loop over `n=1:numel(B)`.
- Line 55: `for` loop over `k=1:numel(B)`.
- Line 86: `for` loop over `n=1:numel(T)`.

### Key state/data transformations

- Lines 20-22: computes `result` using `result=new_test_result('kernel/operator_elementary_suite', 'Elementary operator generators', 'low-level operator constructors must satisfy their defining algebraic ident…`.
- Lines 25: computes `S` using `S=pauli(3)`.
- Lines 27: computes `'angular momentum components obey [Sx,Sy]` using `'angular momentum components obey [Sx,Sy]=i Sz')`.
- Lines 34: computes `W` using `W=weyl(4)`.
- Lines 43: computes `B` using `B=boson_mono(3)`.
- Lines 53: computes `gram` using `gram=zeros(numel(B),numel(B))`.
- Lines 56: computes `gram(n,k)` using `gram(n,k)=trace(full(B{n}'*B{k}))`.
- Lines 63: computes `A` using `A=sin_tran(4)`.
- Lines 64: computes `ref` using `ref=sparse(1,4,1,4,4)`.
- Lines 72: computes `ct_z` using `ct_z=spalloc(4,4,2); ct_z(2,2)=0.5; ct_z(3,3)=-0.5`.
- Lines 73: computes `ct_p` using `ct_p=spalloc(4,4,1); ct_p(2,3)=1`.
- Lines 74: computes `ct_m` using `ct_m=spalloc(4,4,1); ct_m(3,2)=1`.
- Lines 83: computes `T` using `T=irr_sph_ten(3,2)`.
- Lines 84: computes `proj` using `proj=[2 1 0 -1 -2]`.
- Lines 85: computes `L` using `L=pauli(3)`.
- Lines 88: computes `'irreducible spherical tensors obey [Lz,T(k,m)]` using `'irreducible spherical tensors obey [Lz,T(k,m)]=m T(k,m)')`.
- Lines 94: computes `S_pos` using `S_pos=stevens(3,2,+1)`.
- Lines 95: computes `S_neg` using `S_neg=stevens(3,2,-1)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks analytic commutation relations, indexing conventions,
- and small explicit matrices for low-level operator constructors.

## Implementation structure

- Tests elementary operator generators in kernel/operators. Syntax:
- result=test_operator_elementary_suite()
- result -regression test result with explanatory messages
- The test checks analytic commutation relations, indexing conventions,
- and small explicit matrices for low-level operator constructors.
- Announce the test target
- State the operator target of the test
- Check spin-one angular momentum commutation and ladder definitions
- Check finite-truncation Weyl algebra away from the unavoidable edge state
- Check bosonic monomial serpentine indexing
- Check Gram-Schmidt orthogonality without imposing normalisation
- Check single-transition basis indexing from the documented 4x4 map

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `pauli()`, `test_close()`, `speye()`, `weyl()`, `boson_mono()`, `boson_ortho()`, `gram()`, `sin_tran()`, `the()`, `spalloc()`, `ct_z()`, `ct_p()`, `ct_m()`, `centrans()`, `complex()`.
