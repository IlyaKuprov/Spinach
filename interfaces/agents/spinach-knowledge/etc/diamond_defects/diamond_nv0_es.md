# etc/diamond_defects/diamond_nv0_es.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_nv0_es.m`
- Signature: `[sys,inter]=diamond_nv0_es(parameters)`
- Total lines: 108

## Purpose

NV0 excited-state spin system for diamond. Syntax: [sys,inter]=diamond_nv0_es(parameters) Magnetic parameters from Felton et al., Phys. Rev. B 77, 081201 (2008), https://doi.org/10.1103/PhysRevB.77.081201

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters is a structure with the following required fields:
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nitrogen -'14N' or '15N'; 14N hyperfine couplings
- are scaled from 15N, with no NQI included
- because none is reported for this state

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- NV0 excited-state spin system for diamond. Syntax:
- [sys,inter]=diamond_nv0_es(parameters)
- Magnetic parameters from Felton et al., Phys. Rev. B 77,
- 081201 (2008), https://doi.org/10.1103/PhysRevB.77.081201
- parameters is a structure with the following required fields:
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nitrogen -'14N' or '15N'; 14N hyperfine couplings
- are scaled from 15N, with no NQI included
- because none is reported for this state
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `zfs2mat()`, `spin()`, `rotmat_align()`, `mat2ias()`, `isstruct()`, `isfield()`, `ischar()`, `ismember()`, `any()`, `strcmp()`.
