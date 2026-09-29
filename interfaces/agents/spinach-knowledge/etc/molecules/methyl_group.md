# etc/molecules/methyl_group.m

- MATLAB implementation: [etc/molecules/methyl_group.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/molecules/methyl_group.m)

**Call:** `xyz = methyl_group(c_xyz, cc_th, cc_ph, phase)` (for example, `xyz=methyl_group([0 0 0],pi/2,0,0)` with Spinach's `euler2dcm` available).

**Inputs:** `c_xyz` is a real numeric three-element row vector in Angstrom. `cc_th` and `cc_ph` are real numeric scalar polar and azimuthal angles for the C–C bond, in radians; `phase` is a real numeric scalar phase in radians for methyl rotation around that bond. The checks enforce real numeric scalar angles and the stated row-vector shape; the source does not impose angle ranges or an explicit finiteness test.

## Geometry construction

The carbon is placed at the supplied `c_xyz`. Three C–H vectors of length 1.050 Angstrom are generated with the tetrahedral polar angle `acos(1/3)` and azimuths `phase`, `2*pi/3 + phase`, and `4*pi/3 + phase`. The routine first rotates this canonical four-atom geometry to orient the C–C bond using `euler2dcm(0,-cc_th,-cc_ph)`, then translates every coordinate by `c_xyz`.

## Output

`xyz` is a 4-by-1 cell array of Cartesian 1-by-3 coordinate row vectors in Angstrom: carbon first, followed by the three hydrogens. This function returns coordinates only; it does not create a Spinach spin-system structure.

**Source reference:** [Spinach Wiki: methyl_group.m](https://spindynamics.org/wiki/index.php?title=methyl_group.m).
