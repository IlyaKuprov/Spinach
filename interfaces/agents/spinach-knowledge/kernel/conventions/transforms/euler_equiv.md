# kernel/conventions/transforms/euler_equiv.m

**MATLAB source:** [kernel/conventions/transforms/euler_equiv.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/euler_equiv.m)
**Spinach Wiki:** [euler_equiv.m](https://spindynamics.org/wiki/index.php?title=euler_equiv.m)

## Purpose and convention

Tests whether two ZYZ active Euler-angle triples represent rotations whose relative geodesic angle is no greater than a tolerance. Euler triples are not compared component-by-component because they are not unique; the matrices produced by [euler2dcm.m](euler2dcm.md) are compared instead. All angles and the tolerance are in radians.

## Inputs

Signature: `answer=euler_equiv(eulers_a,eulers_b,tol)`

- `eulers_a` and `eulers_b`: numeric, real, finite three-element vectors `[alpha beta gamma]` in the ZYZ active convention. Row and column vectors are accepted; each is reshaped to a row internally.
- `tol`: numeric, real, finite scalar with `tol >= 0`.

## Comparison

Let `D_a=euler2dcm(eulers_a)` and `D_b=euler2dcm(eulers_b)`. The implementation forms `D_rel = D_b * D_a^T` and calculates:

- `s = norm([D_rel(3,2)-D_rel(2,3); D_rel(1,3)-D_rel(3,1); D_rel(2,1)-D_rel(1,2)], 2) / 2`
- `c = (trace(D_rel)-1) / 2`, then clamps `c` to `[-1,1]`.
- `theta = atan2(s,c)`.

The scalar logical output is `answer = (theta <= tol)`. Thus equality at the requested tolerance returns true.
