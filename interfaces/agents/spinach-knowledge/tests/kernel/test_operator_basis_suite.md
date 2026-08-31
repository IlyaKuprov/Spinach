# tests/kernel/test_operator_basis_suite.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/tests/kernel/test_operator_basis_suite.m`
- Signature: `result=test_operator_basis_suite()`
- Total lines: 141

## Purpose

Tests operator-basis construction and expansion helpers. Syntax: result=test_operator_basis_suite()

## Physical / mathematical content

- This file belongs to the `tests` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file also defines local helper function(s): `ist_reconstruct()`, `bm_reconstruct()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 17-18: Announce the test target; implemented by `fprintf('TESTING: Operator-basis construction functions\n')`.
- Lines 20-23: State the operator-basis target of the test; implemented by `result=new_test_result('kernel/operator_basis_suite', 'Operator-basis construction functions', 'operator-basis helpers must satisfy their defining algebra and reconstruc…`.
- Lines 25-26: Irreducible spherical tensors obey [Lz,T(k,m)]=m*T(k,m); implemented by `mult=3; L=pauli(mult); T=irr_sph_ten(mult,2); projections=2:-1:-2`.
- Lines 36-38: Stevens rank-one zero-projection operator is Lz; implemented by `result=test_close(result,'stevens rank-one z',stevens(3,1,0),L.z,1e-14,1e-14, 'the rank-one q=0 Stevens operator is the angular momentum z operator')`.
- Lines 40-41: Weyl boson operators obey their defining number-operator commutators; implemented by `W=weyl(4)`.
- Lines 49-50: Bosonic monomials and orthogonalised monomials must have documented structure; implemented by `B=boson_mono(3)`.
- Lines 65-66: Single-transition operators follow serpentine indexing exactly; implemented by `ST=sin_tran(3)`.
- Lines 73-74: Central-transition operators act only on the central two Zeeman levels; implemented by `CTz=centrans(4,'z')`.
- Lines 81-82: IST and bosonic expansion helpers must reconstruct their source operators; implemented by `A=[1 2i;-3 4]`.
- Lines 103-104: Two-spin IST and sparse preallocation dimensions must match the active basis; implemented by `sys.magnet=0`.

### Control flow inferred from the code

- Line 27: `for` loop over `n=1:numel(T)`.
- Line 57: `for` loop over `n=1:numel(Bo)`.
- Line 58: `for` loop over `k=1:(n-1)`.
- Line 67: `for` loop over `n=1:9`.

### Key state/data transformations

- Lines 21-23: computes `result` using `result=new_test_result('kernel/operator_basis_suite', 'Operator-basis construction functions', 'operator-basis helpers must satisfy their defining algebra and reconstruc…`.
- Lines 26: computes `mult` using `mult=3; L=pauli(mult); T=irr_sph_ten(mult,2); projections=2:-1:-2`.
- Lines 32: computes `Tall` using `Tall=irr_sph_ten(mult)`.
- Lines 38: computes `'the rank-one q` using `'the rank-one q=0 Stevens operator is the angular momentum z operator')`.
- Lines 41: computes `W` using `W=weyl(4)`.
- Lines 50: computes `B` using `B=boson_mono(3)`.
- Lines 55: computes `Bo` using `Bo=boson_ortho(3)`.
- Lines 56: computes `max_overlap` using `max_overlap=0`.
- Lines 66: computes `ST` using `ST=sin_tran(3)`.
- Lines 68: computes `[row,col]` using `[row,col]=lin2kq(3,n,1)`.
- Lines 74: computes `CTz` using `CTz=centrans(4,'z')`.
- Lines 75: computes `CTp` using `CTp=centrans(4,'+')`.
- Lines 82: computes `A` using `A=[1 2i;-3 4]`.
- Lines 83: computes `[states,coeffs]` using `[states,coeffs]=oper2ist(A)`.
- Lines 90: computes `P` using `P=zeros(3); P(2,2)=1`.
- Lines 104: computes `sys.magnet` using `sys.magnet=0`.
- Lines 105: computes `sys.isotopes` using `sys.isotopes={'1H','1H'}`.
- Lines 106: computes `inter.zeeman.scalar` using `inter.zeeman.scalar={0,0}`.

### Local helper functions

- Line 126: `ist_reconstruct()` — `function A=ist_reconstruct(mult,states,coeffs)`. Reconstruct a matrix from bosonic monomial coefficients.
  - Representative operation: `T=irr_sph_ten(mult)`.
  - Representative operation: `A=zeros(mult)`.
- Line 135: `bm_reconstruct()` — `function A=bm_reconstruct(nlevels,states,coeffs)`.
  - Representative operation: `B=boson_mono(nlevels)`.
  - Representative operation: `A=zeros(nlevels)`.

## Outputs

- result -regression test result with explanatory messages
- The test checks single-spin tensor bases, bosonic bases, central
- transitions, Stevens operators, single-transition matrices, expansion
- helpers, and sparse preallocation dimensions.

## Implementation structure

- Tests operator-basis construction and expansion helpers. Syntax:
- result=test_operator_basis_suite()
- result -regression test result with explanatory messages
- The test checks single-spin tensor bases, bosonic bases, central
- transitions, Stevens operators, single-transition matrices, expansion
- helpers, and sparse preallocation dimensions.
- Announce the test target
- State the operator-basis target of the test
- Irreducible spherical tensors obey [Lz,T(k,m)]=m*T(k,m)
- Stevens rank-one zero-projection operator is Lz
- Weyl boson operators obey their defining number-operator commutators
- Bosonic monomials and orthogonalised monomials must have documented structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `new_test_result()`, `pauli()`, `irr_sph_ten()`, `test_close()`, `int2str()`, `projections()`, `comm()`, `stevens()`, `weyl()`, `boson_mono()`, `speye()`, `boson_ortho()`, `hdot()`, `sin_tran()`, `lin2kq()`, `centrans()`.
