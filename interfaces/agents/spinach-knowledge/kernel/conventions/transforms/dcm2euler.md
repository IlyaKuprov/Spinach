# kernel/conventions/transforms/dcm2euler.m

- Signature: `[arg1,arg2,arg3]=dcm2euler(dcm)`

## Purpose

Converts directional cosine matrix into Euler angles, ZYZ active convention (rotating the object rather than the axes). Syntax: [alpha,beta,gamma]=dcm2euler(dcm) OR angles=dcm2euler(dcm)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- An eigenvalue problem is solved or analysed, so the file is extracting spectra, stationary states, avoided crossings, or modal structure from the effective Hamiltonian or superoperator.

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
