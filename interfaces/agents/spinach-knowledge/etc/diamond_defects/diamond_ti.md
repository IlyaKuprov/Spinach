# etc/diamond_defects/diamond_ti.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_ti.m`
- Signature: `[sys,inter]=diamond_ti(parameters)`
- Total lines: 186

## Purpose

Titanium-related defect spin system for diamond. Syntax: [sys,inter]=diamond_ti(parameters) Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `diamond_frame_xyz()`, `diamond_frame_alpha()`, `diamond_frame_xz()`, `round()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'n3' or 'ok1'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .titanium -titanium isotope label, or 'none'
- .n_13c -number of reported 13C hyperfine couplings
- to include, from 0 to 2; applies only to OK1

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Titanium-related defect spin system for diamond. Syntax:
- [sys,inter]=diamond_ti(parameters)
- Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),
- parameters is a structure with the following fields:
- .centre -'n3' or 'ok1'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .titanium -titanium isotope label, or 'none'
- .n_13c -number of reported 13C hyperfine couplings
- to include, from 0 to 2; applies only to OK1
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmpi()`, `spin()`, `lower()`, `diamond_frame_alpha()`, `strcmp()`, `diamond_frame_xz()`, `rotmat_align()`, `isfield()`, `diamond_frame_xyz()`, `xaxis()`, `yaxis()`, `zaxis()`, `frame()`, `cosd()`, `sind()`.
