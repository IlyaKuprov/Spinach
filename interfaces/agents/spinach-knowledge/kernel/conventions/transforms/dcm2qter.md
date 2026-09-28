# kernel/conventions/transforms/dcm2qter.m

- Signature: `q=dcm2qter(dcm)`

## Purpose

Converts a direction cosine matrix in the active convention of euler2dcm.m function into a unit quaternion. Syntax: q=dcm2qter(dcm)

## Physical / mathematical content
Represents a 3D rotation matrix as a unit quaternion in the active convention used by `euler2dcm.m`. Because `q` and `-q` represent the same rotation, the returned quaternion is chosen to have a non-negative scalar component `q.u`.

## Numerical / algorithmic content
Computes four Shepperd pivot invariants from the matrix diagonal and selects the largest for a better-conditioned conversion. The selected component is obtained by a square root; the other three are computed from symmetric or antisymmetric matrix entries. The result is sign-adjusted if `q.u < 0` and then divided by its Euclidean norm.

## Parameters / inputs

- dcm -directional cosine matrix, a 3x3 orthogonal
- matrix with unit determinant

## Outputs

- q -structure with four scalar fields q.u, q.i, q.j,
- q.k giving the four components of the quaternion,
- normalised to q.u greater than or equal to zero
- Note: quaternions double-cover rotations; of the two candi-
- dates q and -q this function returns the one with the
- non-negative scalar part.

## Implementation structure
The main function calls the local `grumble` validator, computes the four candidate pivots, and uses a four-case switch to construct `q.u`, `q.i`, `q.j`, and `q.k`. It then resolves the quaternion sign and normalizes the result. `grumble` requires a real numeric 3×3 matrix and checks orthogonality and unit determinant, each with a tolerance of `1e-6`.
