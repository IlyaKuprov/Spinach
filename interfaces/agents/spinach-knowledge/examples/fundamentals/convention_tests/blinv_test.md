# examples/fundamentals/convention_tests/blinv_test.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/convention_tests/blinv_test.m`
- Signature: `blinv_test()`
- Total lines: 43

## Purpose

Internal consistency of Blicharski invariants and their polarisation relationships with inner products of irre- ducible spherical tensor coefficients, second rank.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 9-10: Random tensors; implemented by `A=rand(3,3); B=rand(3,3)`.
- Lines 12-13: Spherical expansions; implemented by `[~,~,Phi_A]=mat2sphten(A)`.
- Lines 16-17: Blicharski invariants; implemented by `[~,Dsq_A]=blinv(A)`.
- Lines 21-23: Blicharski against Phi products, rank 2 self; implemented by `DA=(Phi_A(3)^2-2*Phi_A(2)*Phi_A(4) +2*Phi_A(1)*Phi_A(5))-(2/3)*Dsq_A`.
- Lines 32-35: Blicharski against Phi product, rank 2 cross; implemented by `D_AB=(Phi_A(1)*Phi_B(5)-Phi_A(2)*Phi_B(4)+ Phi_A(3)*Phi_B(3)-Phi_A(4)*Phi_B(2)+ Phi_A(5)*Phi_B(1))-(2/3)*X_AB`.

### Control flow inferred from the code

- Line 26: conditional branch on `(abs(DA)>10*eps)||(abs(DB)>10*eps)`.
- Line 36: conditional branch on `(abs(D_AB)>10*eps)`.

### Key state/data transformations

- Lines 10: computes `A` using `A=rand(3,3); B=rand(3,3)`.
- Lines 13: computes `[~,~,Phi_A]` using `[~,~,Phi_A]=mat2sphten(A)`.
- Lines 14: computes `[~,~,Phi_B]` using `[~,~,Phi_B]=mat2sphten(B)`.
- Lines 17: computes `[~,Dsq_A]` using `[~,Dsq_A]=blinv(A)`.
- Lines 18: computes `[~,Dsq_B]` using `[~,Dsq_B]=blinv(B)`.
- Lines 19: computes `[~,X_AB]` using `[~,X_AB]=blprod(A,B)`.
- Lines 22-23: computes `DA` using `DA=(Phi_A(3)^2-2*Phi_A(2)*Phi_A(4) +2*Phi_A(1)*Phi_A(5))-(2/3)*Dsq_A`.
- Lines 24-25: computes `DB` using `DB=(Phi_B(3)^2-2*Phi_B(2)*Phi_B(4) +2*Phi_B(1)*Phi_B(5))-(2/3)*Dsq_B`.
- Lines 33-35: computes `D_AB` using `D_AB=(Phi_A(1)*Phi_B(5)-Phi_A(2)*Phi_B(4)+ Phi_A(3)*Phi_B(3)-Phi_A(4)*Phi_B(2)+ Phi_A(5)*Phi_B(1))-(2/3)*X_AB`.

## Implementation structure

- Internal consistency of Blicharski invariants and their
- polarisation relationships with inner products of irre-
- ducible spherical tensor coefficients, second rank.
- Random tensors
- Spherical expansions
- Blicharski invariants
- Blicharski against Phi products, rank 2 self
- Blicharski against Phi product, rank 2 cross

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mat2sphten()`, `blinv()`, `blprod()`, `Phi_A()`, `Phi_B()`.
