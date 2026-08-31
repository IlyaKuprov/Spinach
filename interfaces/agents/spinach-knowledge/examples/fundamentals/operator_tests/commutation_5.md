# examples/fundamentals/operator_tests/commutation_5.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/operator_tests/commutation_5.m`
- Signature: `commutation_5()`
- Total lines: 62

## Purpose

Tests of basic product and commutation relations between bosonic operators.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 8-9: Random dimension; implemented by `nlevels=5+randi(20)`.
- Lines 11-12: Accuracy threshold; implemented by `acc=10*nlevels*eps('double')`.
- Lines 14-15: Weyl operators; implemented by `A=weyl(nlevels)`.
- Lines 17-18: Product test; implemented by `err=norm(A.c*A.a-A.n,'fro')/norm(A.n,'fro')`.
- Lines 25-26: Exact commutator test; implemented by `err_a=norm(comm(A.n,A.c)-A.c,'fro')/norm(A.c,'fro')`.
- Lines 34-35: Boundary commutator test; implemented by `C=comm(A.a,A.c); C(end,end)=1`.
- Lines 43-44: Bosonic monomials; implemented by `B=boson_mono(nlevels)`.
- Lines 46-47: Exact commutator test; implemented by `for n=1:numel(B)`.
- Lines 49-50: Physical indices; implemented by `[k,q]=lin2kq(nlevels,n,1)`.
- Lines 52-53: Deviation norm testing; implemented by `err=norm(comm(A.n,B{n})-(k-q)*B{n},'fro')/norm(B{n},'fro')`.

### Control flow inferred from the code

- Line 19: conditional branch on `err>acc`.
- Line 28: conditional branch on `(err_a>acc)||(err_b>acc)`.
- Line 37: conditional branch on `err>acc`.
- Line 47: `for` loop over `n=1:numel(B)`.
- Line 54: conditional branch on `err>acc`.

### Key state/data transformations

- Lines 9: computes `nlevels` using `nlevels=5+randi(20)`.
- Lines 12: computes `acc` using `acc=10*nlevels*eps('double')`.
- Lines 15: computes `A` using `A=weyl(nlevels)`.
- Lines 18: computes `err` using `err=norm(A.c*A.a-A.n,'fro')/norm(A.n,'fro')`.
- Lines 26: computes `err_a` using `err_a=norm(comm(A.n,A.c)-A.c,'fro')/norm(A.c,'fro')`.
- Lines 27: computes `err_b` using `err_b=norm(comm(A.n,A.a)+A.a,'fro')/norm(A.a,'fro')`.
- Lines 35: computes `C` using `C=comm(A.a,A.c); C(end,end)=1`.
- Lines 44: computes `B` using `B=boson_mono(nlevels)`.
- Lines 50: computes `[k,q]` using `[k,q]=lin2kq(nlevels,n,1)`.

## Implementation structure

- Tests of basic product and commutation relations between
- bosonic operators.
- Random dimension
- Accuracy threshold
- Weyl operators
- Product test
- Exact commutator test
- Boundary commutator test
- Bosonic monomials
- Physical indices
- Deviation norm testing

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `randi()`, `eps()`, `weyl()`, `comm()`, `boson_mono()`, `lin2kq()`.
