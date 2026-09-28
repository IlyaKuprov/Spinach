# kernel/utilities/wigner_3j.m

- Signature: `w=wigner_3j(j1,m1,j2,m2,j3,m3)`

Calculates the Wigner 3j-symbol with `j1`, `j2`, `j3` in the top row and `m1`, `m2`, `m3` in the bottom row. Physically inadmissible indices yield zero.

## Inputs and output

- `j1`, `j2`, `j3`: angular-momentum indices, in top-row order.
- `m1`, `m2`, `m3`: magnetic indices, in bottom-row order.
- `w`: resulting 3j-symbol.

All inputs must be numeric, real, scalar integers or half-integers; otherwise the function raises an error. The result is calculated from a Clebsch–Gordan coefficient as `w = (-1)^(-m3+j1+j2) * clebsch_gordan(j3,-m3,j1,m1,j2,m2) / sqrt(2*j3+1)`.

Contact: ilya.kuprov@weizmann.ac.il

Source: <https://spindynamics.org/wiki/index.php?title=wigner_3j.m>