# examples/fundamentals/convention_tests/rotations_1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/examples/fundamentals/convention_tests/rotations_1.m`
- Signature: `rotations_1()`
- Total lines: 129

## Purpose

Tests the internal consistency of kernel rotation functions.

## Physical / mathematical content

- Fundamentals examples. These are unit tests, convention checks, and pedagogical demonstrations of operator algebra, perturbation theory, tensor conventions, symmetry, quadrature, and numerical differentiation.

## Numerical / algorithmic content


## Code-derived implementation details

### Comment-guided execution stages

- Lines 7-8: Generate a random symmetric traceless 3x3 matrix; implemented by `A=rand(3); A=A+A'; A=A-trace(A)/3`.
- Lines 10-11: Generate a random set of Euler angles; implemented by `eulers=rand(1,3).*[2*pi pi 2*pi]`.
- Lines 13-16: % Test 1: euler2dcm, wigner, mat2sphten; implemented by `[~,~,rank2_a]=mat2sphten(euler2dcm(eulers)*A*euler2dcm(eulers)')`.
- Lines 15-16: DCM rotation followed by a transformation into irreducible components; implemented by `[~,~,rank2_a]=mat2sphten(euler2dcm(eulers)*A*euler2dcm(eulers)')`.
- Lines 18-19: Transformation into irreducible components followed by a Wigner rotation; implemented by `[~,~,rank2]=mat2sphten(A); rank2_b=wigner(2,eulers(1),eulers(2),eulers(3))*rank2`.
- Lines 21-22: Check the difference; implemented by `if norm(rank2_a-rank2_b,2)<1e-10`.
- Lines 28-31: % Test 2: euler2dcm, dcm2euler; implemented by `R=euler2dcm(eulers)`.
- Lines 30-31: Transforms Euler angles into DCM; implemented by `R=euler2dcm(eulers)`.
- Lines 33-34: Transform the DCM back into Euler angles; implemented by `new_eulers=dcm2euler(R)`.
- Lines 36-37: Check the difference; implemented by `if norm(eulers-new_eulers,2)<1e-3`.
- Lines 43-46: % Test 3: euler2dcm, wigner, dcm2wigner; implemented by `W_a=dcm2wigner(euler2dcm(eulers))`.
- Lines 45-46: Transforms Euler angles into DCM, then DCM to Wigner matrix; implemented by `W_a=dcm2wigner(euler2dcm(eulers))`.
- Lines 48-49: Transform Euler angles to Wigner matrix; implemented by `W_b=wigner(2,eulers(1),eulers(2),eulers(3))`.
- Lines 51-52: Check the difference; implemented by `if norm(W_a-W_b,2)<1e-10`.
- Lines 58-61: % Test 4: mat2sphten, sphten2mat; implemented by `[rank0,rank1,rank2]=mat2sphten(A)`.
- Lines 60-61: Transform the matrix into spherical tensors; implemented by `[rank0,rank1,rank2]=mat2sphten(A)`.
- Lines 63-64: Transform back; implemented by `B=sphten2mat(rank0,rank1,rank2)`.
- Lines 66-67: Check the difference; implemented by `if norm(A-B,2)<1e-10`.

### Control flow inferred from the code

- Line 22: conditional branch on `norm(rank2_a-rank2_b,2)<1e-10`.
- Line 37: conditional branch on `norm(eulers-new_eulers,2)<1e-3`.
- Line 52: conditional branch on `norm(W_a-W_b,2)<1e-10`.
- Line 67: conditional branch on `norm(A-B,2)<1e-10`.
- Line 79: conditional branch on `norm(R-qter2dcm(dcm2qter(R)),2)<1e-10`.
- Line 99: conditional branch on `norm([q1.u-q2.u q1.i-q2.i q1.j-q2.j q1.k-q2.k],2)<1e-10`.
- Line 122: conditional branch on `norm(R1-R2,2)<1e-10`.

### Key state/data transformations

- Lines 8: computes `A` using `A=rand(3); A=A+A'; A=A-trace(A)/3`.
- Lines 11: computes `eulers` using `eulers=rand(1,3).*[2*pi pi 2*pi]`.
- Lines 16: computes `[~,~,rank2_a]` using `[~,~,rank2_a]=mat2sphten(euler2dcm(eulers)*A*euler2dcm(eulers)')`.
- Lines 19: computes `[~,~,rank2]` using `[~,~,rank2]=mat2sphten(A); rank2_b=wigner(2,eulers(1),eulers(2),eulers(3))*rank2`.
- Lines 31: computes `R` using `R=euler2dcm(eulers)`.
- Lines 34: computes `new_eulers` using `new_eulers=dcm2euler(R)`.
- Lines 46: computes `W_a` using `W_a=dcm2wigner(euler2dcm(eulers))`.
- Lines 49: computes `W_b` using `W_b=wigner(2,eulers(1),eulers(2),eulers(3))`.
- Lines 61: computes `[rank0,rank1,rank2]` using `[rank0,rank1,rank2]=mat2sphten(A)`.
- Lines 64: computes `B` using `B=sphten2mat(rank0,rank1,rank2)`.
- Lines 88: computes `q1.u` using `q1.u=randn(); q1.i=randn()`.
- Lines 89: computes `q1.j` using `q1.j=randn(); q1.k=randn()`.
- Lines 90: computes `qnorm` using `qnorm=norm([q1.u q1.i q1.j q1.k],2)`.
- Lines 95: computes `[aa_axis,aa_angle]` using `[aa_axis,aa_angle]=qter2anax(q1)`.
- Lines 96: computes `q2` using `q2=anax2qter(aa_axis,aa_angle)`.
- Lines 108: computes `q.u` using `q.u=randn(); q.i=randn()`.
- Lines 109: computes `q.j` using `q.j=randn(); q.k=randn()`.
- Lines 115: computes `R1` using `R1=qter2dcm(q)`.

## Implementation structure

- Tests the internal consistency of kernel rotation functions.
- Generate a random symmetric traceless 3x3 matrix
- Generate a random set of Euler angles
- % Test 1: euler2dcm, wigner, mat2sphten
- DCM rotation followed by a transformation into irreducible components
- Transformation into irreducible components followed by a Wigner rotation
- Check the difference
- % Test 2: euler2dcm, dcm2euler
- Transforms Euler angles into DCM
- Transform the DCM back into Euler angles
- % Test 3: euler2dcm, wigner, dcm2wigner
- Transforms Euler angles into DCM, then DCM to Wigner matrix

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `mat2sphten()`, `euler2dcm()`, `wigner()`, `eulers()`, `dcm2euler()`, `dcm2wigner()`, `sphten2mat()`, `qter2dcm()`, `dcm2qter()`, `qter2anax()`, `anax2qter()`, `anax2dcm()`.
