# etc/diamond_defects/diamond_ni.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_ni.m`
- Signature: `[sys,inter]=diamond_ni(parameters)`
- Total lines: 232

## Purpose

Nickel-related defect spin system for diamond. Syntax: [sys,inter]=diamond_ni(parameters) W8 magnetic parameters from Ludwig and Woodbury, Phys. Rev. B 41, 3905 (1990), https://doi.org/10.1103/PhysRevB.41.3905 The W8 quartet entry assumes zero ZFS: no ZFS parameters were reported in the cited data, and off-central transitions are treated as unresolved rather than explicitly modelled. Other nickel-centre table values 

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `diamond_frame_xyz()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'w8', 'ne1', 'ne2', 'ne3', 'ne4', 'ne5',
- 'ne8', 'ab1', 'ab2', 'ab3', 'ab4', 'ab5',
- 'nol1', or 'nirim5'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nickel -'61Ni', 'none', or another nickel isotope;
- required when .centre is 'w8'
- .n_13c -number of reported 13C hyperfine couplings
- to include, from 0 to 4; applies only to W8

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Nickel-related defect spin system for diamond. Syntax:
- [sys,inter]=diamond_ni(parameters)
- W8 magnetic parameters from Ludwig and Woodbury, Phys. Rev. B 41,
- 3905 (1990), https://doi.org/10.1103/PhysRevB.41.3905
- The W8 quartet entry assumes zero ZFS: no ZFS parameters were
- reported in the cited data, and off-central transitions are treated
- as unresolved rather than explicitly modelled.
- Other nickel-centre table values from Nadolinny et al., Crystals
- 7, 237 (2017), https://doi.org/10.3390/cryst7080237
- parameters is a structure with the following fields:
- .centre -'w8', 'ne1', 'ne2', 'ne3', 'ne4', 'ne5',
- 'ne8', 'ab1', 'ab2', 'ab3', 'ab4', 'ab5',

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmpi()`, `spin()`, `lower()`, `strcmp()`, `nuclei()`, `cosd()`, `sind()`, `cross()`, `diamond_frame_xyz()`, `avalues()`, `zfs2mat()`, `rotmat_align()`, `mat2ias()`, `isfield()`, `xaxis()`.
