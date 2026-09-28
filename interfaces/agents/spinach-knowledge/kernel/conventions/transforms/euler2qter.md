# kernel/conventions/transforms/euler2qter.m

- Signature: `q=euler2qter(arg1,arg2,arg3)`

## Purpose

Converts ZYZ active Euler angles to a unit quaternion in the active convention, matching the rotation represented by `euler2dcm`. Call with three angles, `q=euler2qter(alpha,beta,gamma)`, or one three-element vector, `q=euler2qter([alpha beta gamma])`.

## Parameters / inputs

- `alpha,beta,gamma`: Euler angles in radians; each may be a real scalar or an equal-length column vector.

## Output

- `q`: structure with fields `u`, `i`, `j`, and `k`, containing the four quaternion components. For vector inputs, each field is a column vector. The quaternion converts to the corresponding rotation matrix with `qter2dcm`.

## Quaternion components

`q.u = cos(beta/2) cos((alpha+gamma)/2)`

`q.i = sin(beta/2) sin((gamma-alpha)/2)`

`q.j = sin(beta/2) cos((gamma-alpha)/2)`

`q.k = cos(beta/2) sin((alpha+gamma)/2)`

Source: [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=euler2qter.m)
