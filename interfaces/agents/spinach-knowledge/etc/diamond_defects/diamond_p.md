# etc/diamond_defects/diamond_p.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_p.m`
- Signature: `[sys,inter]=diamond_p(parameters)`
- Total lines: 170

## Purpose

Phosphorus-related defect spin system for diamond. Syntax: [sys,inter]=diamond_p(parameters) Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'ma1', 'np1', 'np2', 'np3', 'np4', 'np5',
- 'np6', 'np8', or 'np9'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .include_13c -include the reported 13C hyperfine coupling;
- applies only to MA1 and defaults to false

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Phosphorus-related defect spin system for diamond. Syntax:
- [sys,inter]=diamond_p(parameters)
- Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),
- parameters is a structure with the following fields:
- .centre -'ma1', 'np1', 'np2', 'np3', 'np4', 'np5',
- 'np6', 'np8', or 'np9'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .include_13c -include the reported 13C hyperfine coupling;
- applies only to MA1 and defaults to false
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isstruct()`, `isfield()`, `grumble()`, `spin()`, `lower()`, `rotmat_align()`, `ischar()`, `ismember()`, `islogical()`, `isscalar()`, `strcmpi()`.
