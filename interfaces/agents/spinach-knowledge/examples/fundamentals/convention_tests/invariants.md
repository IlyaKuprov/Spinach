# examples/fundamentals/convention_tests/invariants.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/convention_tests/invariants.m`
- Signature: `invariants()`
- Total lines: 60

## Purpose

A test of Equation 3 in http://dx.doi.org/10.1002/chem.200902300

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 7-8: Rhombic test; implemented by `for n=1:100`.
- Lines 10-11: Random symmetric matrix; implemented by `A=randn(3); A=A+A'`.
- Lines 13-14: Eigenvalues; implemented by `[~,D]=eig(A); D=diag(D); D=D(randperm(3))`.
- Lines 16-17: Axiality and rhombicity; implemented by `Ax=2*D(3)-(D(1)+D(2)); Rh=D(1)-D(2)`.
- Lines 19-20: Standard relaxation theory invariant; implemented by `P=(Ax^2+3*Rh^2)/6`.
- Lines 22-23: Neat eigenvalue form; implemented by `Q=(2/3)*(D(1)^2+D(2)^2+D(3)^2-D(1)*D(2)-D(1)*D(3)-D(2)*D(3))`.
- Lines 25-26: Difference and diagnostics; implemented by `if abs(P-Q)>1e-6, error('Rhombic test failed.'); end`.
- Lines 30-31: Successful completion message; implemented by `if n==100, disp('Rhombic test passed.'); end`.
- Lines 33-34: Axial test; implemented by `for n=1:100`.
- Lines 36-37: Random axial eigenvalues; implemented by `D=randn(2,1); D(3)=D(2)`.
- Lines 39-40: Zero trace; implemented by `D=D-sum(D)/3`.
- Lines 42-43: Standard invarinat; implemented by `P=(D(1)-D(2))^2`.
- Lines 45-46: Random shuffle; implemented by `D=D(randperm(3))`.
- Lines 48-49: Neat invariant; implemented by `Q=D(1)^2+D(2)^2+D(3)^2-D(1)*D(2)-D(1)*D(3)-D(2)*D(3)`.
- Lines 51-52: Difference and diagnostics; implemented by `if abs(P-Q)>1e-6, error('Axial test failed.'); end`.
- Lines 56-57: Successful completion message; implemented by `if n==100, disp('Axial test passed.'); end`.

### Control flow inferred from the code

- Line 8: `for` loop over `n=1:100`.
- Line 26: conditional branch on `abs(P-Q)>1e-6, error('Rhombic test failed.'); end`.
- Line 31: conditional branch on `n==100, disp('Rhombic test passed.'); end`.
- Line 34: `for` loop over `n=1:100`.
- Line 52: conditional branch on `abs(P-Q)>1e-6, error('Axial test failed.'); end`.
- Line 57: conditional branch on `n==100, disp('Axial test passed.'); end`.

### Key state/data transformations

- Lines 11: computes `A` using `A=randn(3); A=A+A'`.
- Lines 14: computes `[~,D]` using `[~,D]=eig(A); D=diag(D); D=D(randperm(3))`.
- Lines 17: computes `Ax` using `Ax=2*D(3)-(D(1)+D(2)); Rh=D(1)-D(2)`.
- Lines 20: computes `P` using `P=(Ax^2+3*Rh^2)/6`.
- Lines 23: computes `Q` using `Q=(2/3)*(D(1)^2+D(2)^2+D(3)^2-D(1)*D(2)-D(1)*D(3)-D(2)*D(3))`.
- Lines 37: computes `D` using `D=randn(2,1); D(3)=D(2)`.

## Implementation structure

- A test of Equation 3 in http://dx.doi.org/10.1002/chem.200902300
- Rhombic test
- Random symmetric matrix
- Eigenvalues
- Axiality and rhombicity
- Standard relaxation theory invariant
- Neat eigenvalue form
- Difference and diagnostics
- Successful completion message
- Axial test
- Random axial eigenvalues
- Zero trace

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `randperm()`.
