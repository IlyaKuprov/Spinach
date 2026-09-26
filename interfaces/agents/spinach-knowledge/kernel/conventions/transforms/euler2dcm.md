# kernel/conventions/transforms/euler2dcm.m

- Signature: `R=euler2dcm(arg1,arg2,arg3)`

## Purpose

Converts Euler angles (ZYZ active convention) into a direction cosine matrix. Syntax: R=euler2dcm(alpha,beta,gamma) OR R=euler2dcm([alpha beta gamma])

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- alpha,beta,gamma -Euler angles in radians (ZYZ
- active convention)

## Outputs

- R -direction cosine matrix
- Note: the resulting rotation matrix is to be used as follows:
- v=R*v (for 3x1 vectors)
- A=R*A*R' (for 3x3 interaction tensors)

## Implementation structure

- Converts Euler angles (ZYZ active convention) into a direction
- cosine matrix. Syntax:
- R=euler2dcm(alpha,beta,gamma)
- R=euler2dcm([alpha beta gamma])
- alpha,beta,gamma - Euler angles in radians (ZYZ
- active convention)
- R - direction cosine matrix
- Note: the resulting rotation matrix is to be used as follows:
- v=R*v (for 3x1 vectors)
- A=R*A*R' (for 3x3 interaction tensors)
- Adapt to the input style
- Assume that a single input is a 3-vector
