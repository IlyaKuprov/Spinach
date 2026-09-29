# kernel/conventions/transforms/euler2qter.m

**MATLAB source:** [kernel/conventions/transforms/euler2qter.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/euler2qter.m)
**Spinach Wiki:** [euler2qter.m](https://spindynamics.org/wiki/index.php?title=euler2qter.m)

## Purpose and convention

Converts ZYZ active Euler angles, in radians, to the unit quaternion representing the same active rotation as [euler2dcm.m](euler2dcm.md). The quaternion can be converted to that direction-cosine matrix with [qter2dcm.m](qter2dcm.md).

## Inputs

Signature: `q=euler2qter(arg1,arg2,arg3)`

The function accepts either one argument, interpreted by indexing its first three elements as `[alpha beta gamma]`, or three separate arguments `alpha`, `beta`, and `gamma`. In the three-argument form each value must be numeric, real, and a column (a scalar is also accepted); the three inputs must have equal element counts. The one-argument branch does not validate vector shape or length before indexing: fewer than three elements fail during indexing, and elements after the first three are unused. No finite-value check is applied to the angles.

Angles are in radians and use the ZYZ active convention.

## Output and conversion

Returns a structure with fields `q.u`, `q.i`, `q.j`, and `q.k`. For vector-valued three-argument inputs, each field is a column vector; scalar inputs produce scalar fields.

The source computes:

- `q.u = cos(beta/2) * cos((alpha+gamma)/2)`
- `q.i = sin(beta/2) * sin((gamma-alpha)/2)`
- `q.j = sin(beta/2) * cos((gamma-alpha)/2)`
- `q.k = cos(beta/2) * sin((alpha+gamma)/2)`

The returned quaternion is in the active convention and corresponds to the rotation represented by `euler2dcm(alpha,beta,gamma)`.
