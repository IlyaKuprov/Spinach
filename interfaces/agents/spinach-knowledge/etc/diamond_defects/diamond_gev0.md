# etc/diamond_defects/diamond_gev0.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_gev0.m`
- Signature: `[sys,inter]=diamond_gev0(parameters)`
- Total lines: 110

## Purpose

GeV0 spin system for diamond. Syntax: [sys,inter]=diamond_gev0(parameters) Magnetic parameters from Nadolinny et al., Phys. Status Solidi A 213, 2623 (2016), https://doi.org/10.1002/pssa.201600211

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters is a structure with the following fields:
- .germanium -'73Ge', 'none', or another germanium isotope
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- GeV0 spin system for diamond. Syntax:
- [sys,inter]=diamond_gev0(parameters)
- Magnetic parameters from Nadolinny et al., Phys. Status Solidi A
- 213, 2623 (2016), https://doi.org/10.1002/pssa.201600211
- parameters is a structure with the following fields:
- .germanium -'73Ge', 'none', or another germanium isotope
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
- Check input count
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `zfs2mat()`, `strcmp()`, `rotmat_align()`, `mat2ias()`, `isfield()`, `isstruct()`, `ischar()`, `ismember()`.
