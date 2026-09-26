# etc/molecules/methyl_group.m

- Signature: `xyz=methyl_group(c_xyz,cc_th,cc_ph,phase)`

## Purpose

Coordinates for the four atoms of a methyl group. Syntax: xyz=methyl_group(c_xyz,cc_th,cc_ph,phase)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- c_xyz -coordinates of C, row vector, Angstrom
- cc_th -polar theta angle of the C-C bond, radians
- cc_ph -polar phi angle of the C-C bond, radians
- phase -phase of the methyl group with respect to
- its rotation around the C-C bond, radians

## Outputs

- xyz -a column cell array of Cartesian XYZ row
- vectors; carbon is the first atom

## Implementation structure

- Coordinates for the four atoms of a methyl group. Syntax:
- xyz=methyl_group(c_xyz,cc_th,cc_ph,phase)
- c_xyz -coordinates of C, row vector, Angstrom
- cc_th -polar theta angle of the C-C bond, radians
- cc_ph -polar phi angle of the C-C bond, radians
- phase -phase of the methyl group with respect to
- its rotation around the C-C bond, radians
- xyz -a column cell array of Cartesian XYZ row
- vectors; carbon is the first atom
- Check consistency
- Generate a canonical methyl group
- Rotate the CC bond
