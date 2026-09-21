# etc/diamond_defects/diamond_p1.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_p1.m`
- Signature: `[sys,inter]=diamond_p1(parameters)`
- Total lines: 116

## Purpose

P1 centre spin system for diamond. Syntax: [sys,inter]=diamond_p1(parameters) Magnetic parameters from: Nir-Arad et al. Phys. Chem. Chem. Phys. 26, 27633 (2024), <https://doi.org/10.1039/d4cp03055a>, and Smith et al. Phys. Rev. 115, 1546 (1959),

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- a structure (parameters.*) with the following fields:
- .orientation -'111', '110', or '100' crystal
- plane normal aligned with the
- magnetic field, def. is '111'
- .nitrogen -'14N' or '15N', default is '14N'

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- P1 centre spin system for diamond. Syntax:
- [sys,inter]=diamond_p1(parameters)
- Magnetic parameters from: Nir-Arad et al. Phys. Chem. Chem. Phys. 26,
- 27633 (2024), <https://doi.org/10.1039/d4cp03055a>, and
- Smith et al. Phys. Rev. 115, 1546 (1959),
- a structure (parameters.*) with the following fields:
- .orientation -'111', '110', or '100' crystal
- plane normal aligned with the
- magnetic field, def. is '111'
- .nitrogen -'14N' or '15N', default is '14N'
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `rotmat_align()`, `zfs2mat()`, `isstruct()`, `ischar()`.
