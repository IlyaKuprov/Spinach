# etc/diamond_defects/diamond_co.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_co.m`
- Signature: `[sys,inter]=diamond_co(parameters)`
- Total lines: 125

## Purpose

Cobalt-related defect spin system for diamond. Syntax: [sys,inter]=diamond_co(parameters) Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'o4' or 'nlo2'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Cobalt-related defect spin system for diamond. Syntax:
- [sys,inter]=diamond_co(parameters)
- Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),
- parameters is a structure with the following fields:
- .centre -'o4' or 'nlo2'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
- Check input count
- Check consistency
- Set field-unit conversion constants

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `lower()`, `cross()`, `cosd()`, `sind()`, `xaxis()`, `yaxis()`, `zaxis()`, `frame()`, `rotmat_align()`, `isstruct()`, `isfield()`, `ischar()`, `any()`, `strcmpi()`.
