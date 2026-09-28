# kernel/conventions/transforms/dcm2euler.m

- Signature: `[arg1,arg2,arg3]=dcm2euler(dcm)`

## Purpose

Converts directional cosine matrix into Euler angles, ZYZ active convention (rotating the object rather than the axes). Syntax: [alpha,beta,gamma]=dcm2euler(dcm) OR angles=dcm2euler(dcm)

## Physical / mathematical content
The input is a 3×3 directional cosine matrix representing a rotation. The output uses the ZYZ active Euler-angle convention, which rotates the object rather than the axes. For an imperfect input matrix, the function obtains the nearest proper rotation in the Frobenius norm before extracting its angles.

## Numerical / algorithmic content
The function builds a 4×4 Davenport matrix from the input DCM and takes the eigenvector associated with its largest eigenvalue as the quaternion of the nearest proper rotation. It converts that quaternion to Euler angles with `qter2euler`, wraps alpha and gamma using `mod(angle,2*pi)`, and checks the result by reconstructing a DCM. If the reconstruction differs from the input by more than `1e-2` in the matrix 1-norm, it displays both matrices and raises an error.

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
The function first validates that `dcm` is a real, finite 3×3 numeric matrix. Orthogonality and determinant deviations above `1e-6` produce warnings; deviations above `1e-2` produce errors. It then constructs the Davenport matrix, selects its dominant eigenvector, converts the resulting quaternion to ZYZ active Euler angles, wraps alpha and gamma, and checks the reconstructed DCM. With one or no requested outputs it returns `[alpha beta gamma]` as a row vector; with three outputs it returns the angles separately. Any other number of requested outputs raises an error.
