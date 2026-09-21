# etc/diamond_defects/diamond_r2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_r2.m`
- Signature: `[sys,inter]=diamond_r2(parameters)`
- Total lines: 91

## Purpose

R2 self-interstitial spin system for diamond. Syntax: [sys,inter]=diamond_r2(parameters) Magnetic parameters from Hunt et al., Phys. Rev. B 61, 3863 (2000), https://doi.org/10.1103/PhysRevB.61.3863

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters is a structure with the following fields:
- .d_sign -sign of D
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- R2 self-interstitial spin system for diamond. Syntax:
- [sys,inter]=diamond_r2(parameters)
- Magnetic parameters from Hunt et al., Phys. Rev. B 61,
- 3863 (2000), https://doi.org/10.1103/PhysRevB.61.3863
- parameters is a structure with the following fields:
- .d_sign -sign of D
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
- Check input count
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `zfs2mat()`, `rotmat_align()`, `mat2ias()`, `isstruct()`, `isfield()`, `ischar()`, `any()`, `strcmp()`, `isscalar()`.
