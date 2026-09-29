# kernel/utilities/dihedral.m

## Purpose

Computes the dihedral angle between vectors specified by four sets of atomic coordinates, for atoms assumed to be bonded as A-B-C-D.

Source: <https://github.com/IlyaKuprov/Spinach/blob/main/kernel/utilities/dihedral.m>

## Behaviour

- Syntax: `phi=dihedral(A,B,C,D)`.
- The function first validates its arguments via an internal consistency check (`grumble`), which errors with the message `'the arguments must be 3-element row vectors of real numbers.'` if any argument is non-numeric, non-real, does not contain exactly 3 elements, or is not a row vector.
- Unit direction vectors are formed along the three bonds: `b1=(B-A)/norm(B-A,2)`, `b2=(C-B)/norm(C-B,2)`, `b3=(D-C)/norm(D-C,2)`.
- The dihedral angle is computed as `phi=180*atan2(dot(norm(b2,2)*b1,cross(b2,b3)),dot(cross(b1,b2),cross(b2,b3)))/pi`, i.e. via a two-argument arctangent of a dot product against a cross product, converted from radians to degrees by the factor `180/pi`.

## Inputs and outputs

Inputs:

- `A` — row vector of Cartesian coordinates for atom A.
- `B` — row vector of Cartesian coordinates for atom B.
- `C` — row vector of Cartesian coordinates for atom C.
- `D` — row vector of Cartesian coordinates for atom D.

Each must be a 3-element row vector of real numbers.

Output:

- `phi` — dihedral angle, in degrees.

## References

- Spinach Wiki: <https://spindynamics.org/wiki/index.php?title=dihedral.m>
