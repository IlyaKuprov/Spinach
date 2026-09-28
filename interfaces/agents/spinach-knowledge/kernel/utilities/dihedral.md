# kernel/utilities/dihedral.m

- Signature: `phi=dihedral(A,B,C,D)`

## Purpose

Computes the dihedral angle for four atoms assumed to be bonded in the order A-B-C-D.

## Physical / mathematical content

The result is a dihedral angle in degrees.

## Numerical / algorithmic content

The function normalizes the three successive bond vectors `B-A`, `C-B`, and `D-C`, then uses their dot and cross products in `atan2` to calculate the angle. It converts the result from radians to degrees.

## Parameters / inputs

- A -row vector of cartesian coordinates
- for atom A
- B -row vector of cartesian coordinates
- for atom B
- C -row vector of cartesian coordinates
- for atom C
- D -row vector of cartesian coordinates
- for atom D

## Outputs

- phi -dihedral angle, degrees

## Implementation structure

The function first checks that each argument is a real, numeric, three-element row vector, then computes the dihedral angle.

Source reference: <https://spindynamics.org/wiki/index.php?title=dihedral.m>