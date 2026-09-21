# etc/diamond_defects/diamond_n_inter.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_n_inter.m`
- Signature: `[sys,inter]=diamond_n_inter(parameters)`
- Total lines: 131

## Purpose

Nitrogen interstitial spin system for diamond. Syntax: [sys,inter]=diamond_n_inter(parameters) Magnetic parameters from Felton et al., J. Phys. Condens. Matter 21, 364212 (2009), https://doi.org/10.1088/0953-8984/21/36/364212

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `diamond_frame_xyz()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'war9' or 'war10'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nitrogen -'14N' or '15N'

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Nitrogen interstitial spin system for diamond. Syntax:
- [sys,inter]=diamond_n_inter(parameters)
- Magnetic parameters from Felton et al., J. Phys. Condens. Matter
- 21, 364212 (2009), https://doi.org/10.1088/0953-8984/21/36/364212
- parameters is a structure with the following fields:
- .centre -'war9' or 'war10'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nitrogen -'14N' or '15N'
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
- Check input count

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `diamond_frame_xyz()`, `sind()`, `cosd()`, `lower()`, `spin()`, `rotmat_align()`, `xaxis()`, `yaxis()`, `zaxis()`, `frame()`, `isstruct()`, `isfield()`, `ischar()`, `any()`, `strcmp()`.
