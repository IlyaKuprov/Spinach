# examples/fundamentals/derivative_tests/auxmat_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/derivative_tests/auxmat_test.m`
- Signature: `auxmat_test()`
- Total lines: 61

## Purpose

Testing IK's favourite equation numerically -auxiliary matrix expression against a high-accuracy finite diffe- rence approximation.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- Time propagation is explicit. In Spinach this usually means repeated application of matrix exponentials or propagator factorizations to density operators or state vectors in Hilbert/Liouville/Fokker-Planck space.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Random matrices, arbitrary function; implemented by `dim=randi(20); f=@(A)expm(A)`.
- Lines 15-16: Some weird ass function and its derivatives; implemented by `L=@(a,b,c)(a*A+c*cos(b)*B+(a*c^4)*C+a*b*c*B*C)`.
- Lines 21-22: Real coefficients and fin. diff. increment; implemented by `x=randn(1); y=randn(1); z=randn(1); h=sqrt(eps)`.
- Lines 24-25: First parameter; implemented by `P=f([ L(x,y,z) dL_da(x,y,z)`.
- Lines 36-37: Second parameter; implemented by `P=f([ L(x,y,z) dL_db(x,y,z)`.
- Lines 48-49: Third parameter; implemented by `P=f([ L(x,y,z) dL_dc(x,y,z)`.

### Control flow inferred from the code

- Line 30: conditional branch on `norm(P-Q,2)/norm(P,2)>1e-3`.
- Line 42: conditional branch on `norm(P-Q,2)/norm(P,2)>1e-3`.
- Line 54: conditional branch on `norm(P-Q,2)/norm(P,2)>1e-3`.

### Key state/data transformations

- Lines 10: computes `dim` using `dim=randi(20); f=@(A)expm(A)`.
- Lines 11: computes `A` using `A=randn(dim)+1i*randn(dim); A=A/norm(A,2)`.
- Lines 12: computes `B` using `B=randn(dim)+1i*randn(dim); B=B/norm(B,2)`.
- Lines 13: computes `C` using `C=randn(dim)+1i*randn(dim); C=C/norm(B,2)`.
- Lines 16: computes `L` using `L=@(a,b,c)(a*A+c*cos(b)*B+(a*c^4)*C+a*b*c*B*C)`.
- Lines 17: computes `dL_da` using `dL_da=@(a,b,c)(A+(c^4)*C+b*c*B*C)`.
- Lines 18: computes `dL_db` using `dL_db=@(a,b,c)(-c*sin(b)*B+a*c*B*C)`.
- Lines 19: computes `dL_dc` using `dL_dc=@(a,b,c)(cos(b)*B+(4*a*c^3)*C+a*b*B*C)`.
- Lines 22: computes `x` using `x=randn(1); y=randn(1); z=randn(1); h=sqrt(eps)`.
- Lines 25: computes `P` using `P=f([ L(x,y,z) dL_da(x,y,z)`.
- Lines 28-29: computes `Q` using `Q=2*(f(L(x+h,y,z))-f(L(x-h,y,z)))/(3*h)- (f(L(x+2*h,y,z))-f(L(x-2*h,y,z)))/(12*h)`.

## Implementation structure

- Testing IK's favourite equation numerically -auxiliary
- matrix expression against a high-accuracy finite diffe-
- rence approximation.
- Random matrices, arbitrary function
- Some weird ass function and its derivatives
- Real coefficients and fin. diff. increment
- First parameter
- Second parameter
- Third parameter

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `randi()`, `dL_da()`, `dL_db()`, `dL_dc()`.
