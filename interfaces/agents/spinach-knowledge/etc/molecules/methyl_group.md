# etc/molecules/methyl_group.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/molecules/methyl_group.m`
- Signature: `xyz=methyl_group(c_xyz,cc_th,cc_ph,phase)`
- Total lines: 70

## Purpose

Coordinates for the four atoms of a methyl group. Syntax: xyz=methyl_group(c_xyz,cc_th,cc_ph,phase)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `acos()`, `euler2dcm()`, `mat2cell()`, `isvector()`, `isrow()`, `isscalar()`.
