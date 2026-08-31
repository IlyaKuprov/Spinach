# kernel/conventions/transforms/dcm2euler.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/dcm2euler.m`
- Signature: `[arg1,arg2,arg3]=dcm2euler(dcm)`
- Total lines: 114

## Purpose

Converts directional cosine matrix into Euler angles, ZYZ active convention (rotating the object rather than the axes). Syntax: [alpha,beta,gamma]=dcm2euler(dcm) OR angles=dcm2euler(dcm)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 41-42: Check consistency; implemented by `grumble(dcm)`.
- Lines 44-45: Build the Davenport matrix; implemented by `K=[dcm(1,1)+dcm(2,2)+dcm(3,3) dcm(3,2)-dcm(2,3) dcm(1,3)-dcm(3,1) dcm(2,1)-dcm(1,2)`.
- Lines 50-51: Get the quaternion of the nearest rotation; implemented by `[evecs,evals]=eig(K,'vector'); [~,best]=max(evals)`.
- Lines 55-56: Extract the Euler angles; implemented by `[alpha,beta,gamma]=qter2euler(q)`.
- Lines 58-59: Wrap alpha and gamma into [0,2*pi]; implemented by `alpha=mod(alpha,2*pi); gamma=mod(gamma,2*pi)`.
- Lines 61-62: Make sure the result is good enough and bomb out if not; implemented by `if norm(dcm-euler2dcm(alpha,beta,gamma),1)>1e-3`.
- Lines 67-68: Adapt to the output style; implemented by `if nargout==1||nargout==0`.

### Control flow inferred from the code

- Line 62: conditional branch on `norm(dcm-euler2dcm(alpha,beta,gamma),1)>1e-3`.
- Line 68: conditional branch on `nargout==1||nargout==0`.

### Key state/data transformations

- Lines 45: computes `K` using `K=[dcm(1,1)+dcm(2,2)+dcm(3,3) dcm(3,2)-dcm(2,3) dcm(1,3)-dcm(3,1) dcm(2,1)-dcm(1,2)`.
- Lines 51: computes `[evecs,evals]` using `[evecs,evals]=eig(K,'vector'); [~,best]=max(evals)`.
- Lines 52: computes `q.u` using `q.u=evecs(1,best); q.i=evecs(2,best)`.
- Lines 53: computes `q.j` using `q.j=evecs(3,best); q.k=evecs(4,best)`.
- Lines 56: computes `[alpha,beta,gamma]` using `[alpha,beta,gamma]=qter2euler(q)`.
- Lines 59: computes `alpha` using `alpha=mod(alpha,2*pi); gamma=mod(gamma,2*pi)`.
- Lines 69: computes `arg1` using `arg1=[alpha beta gamma]`.

### Local helper functions

- Line 79: `grumble()` — `function grumble(dcm)`.
  - Representative operation: `if (~isnumeric(dcm))||(~isreal(dcm))||(~all(size(dcm)==[3 3]))`.
  - Representative operation: `error('DCM must be a real 3x3 matrix.')`.

## Parameters / inputs

- dcm -directional cosine matrix

## Outputs

- alpha, beta, gamma -Euler angles in ZYZ active con-
- vention, radians
- angles -a row vector of Euler angles in
- ZYZ active convention, ordered
- as alpha, beta, gamma, in radians
- Note: the problem of recovering Euler angles from a DCM is, in
- general, ill-posed. This function is a product of consi-
- derable work, it has passed rigorous testing: it either
- returns a correct answer or gives an informative error.
- Note: the angles returned are those of the proper rotation that
- is nearest to the input in the Frobenius norm; that rota-
- tion is found in closed form through the dominant eigen-
- vector of the Davenport matrix (I.Y. Bar-Itzhack, J. Gui-
- dance Control Dyn. 23 (2000) 1085), and the angles are
- extracted from the corresponding quaternion.

## Implementation structure

- Converts directional cosine matrix into Euler angles, ZYZ active
- convention (rotating the object rather than the axes). Syntax:
- [alpha,beta,gamma]=dcm2euler(dcm)
- angles=dcm2euler(dcm)
- dcm -directional cosine matrix
- alpha, beta, gamma -Euler angles in ZYZ active con-
- vention, radians
- angles -a row vector of Euler angles in
- ZYZ active convention, ordered
- as alpha, beta, gamma, in radians
- Note: the problem of recovering Euler angles from a DCM is, in
- general, ill-posed. This function is a product of consi-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `dcm()`, `evecs()`, `qter2euler()`, `euler2dcm()`, `all()`, `any()`.
