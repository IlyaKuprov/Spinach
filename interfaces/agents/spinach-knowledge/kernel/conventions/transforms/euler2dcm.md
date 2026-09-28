# kernel/conventions/transforms/euler2dcm.m

- Signature: `R=euler2dcm(arg1,arg2,arg3)`

## Purpose

Converts Euler angles (ZYZ active convention) into a direction cosine matrix. Syntax: R=euler2dcm(alpha,beta,gamma) OR R=euler2dcm([alpha beta gamma])

## Physical / mathematical content
Uses the active ZYZ Euler-angle convention. The direction cosine matrix is the product R = Rz(alpha) Ry(beta) Rz(gamma), where each factor is a counterclockwise rotation about its indicated axis. Apply R to a 3×1 vector as R*v and to a 3×3 interaction tensor as R*A*R'.

## Numerical / algorithmic content
Constructs the three rotation matrices from the angles' sine and cosine values, then multiplies them in the order R_alpha*R_beta*R_gamma. Inputs must each be numeric, real, and single-element.

## Parameters / inputs
- `alpha`, `beta`, `gamma`: Euler angles in radians, using the active ZYZ convention. Supply them as three separate inputs or together in one three-element vector.

## Outputs

- R -direction cosine matrix
- Note: the resulting rotation matrix is to be used as follows:
- v=R*v (for 3x1 vectors)
- A=R*A*R' (for 3x3 interaction tensors)

## Implementation structure
With one input, the function reads its first three elements as `alpha`, `beta`, and `gamma`; with three inputs, it uses them directly. Other input counts produce an error. It then validates the three angles, constructs the Z, Y, and Z rotation matrices, and returns their product.
