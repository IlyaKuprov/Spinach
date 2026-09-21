# etc/diamond_defects/diamond_nvm_gs.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_nvm_gs.m`
- Signature: `[sys,inter]=diamond_nvm_gs(parameters)`
- Total lines: 127

## Purpose

NV centre ground state spin system for diamond. Syntax: [sys,inter]=diamond_nvm_gs(parameters) Magnetic parameters from: S. Felton et al., Phys. Rev. B 79, 075203 (2009),

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- the following is needed in the parameters.* structure:
- .orientation -'111', '110', or '100' crystal
- plane normal aligned with the
- magnetic field, def. is '111'
- .nitrogen -'14N' or '15N', default is '14N'

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- NV centre ground state spin system for diamond. Syntax:
- [sys,inter]=diamond_nvm_gs(parameters)
- Magnetic parameters from:
- S. Felton et al., Phys. Rev. B 79, 075203 (2009),
- the following is needed in the parameters.* structure:
- .orientation -'111', '110', or '100' crystal
- plane normal aligned with the
- magnetic field, def. is '111'
- .nitrogen -'14N' or '15N', default is '14N'
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `isfield()`, `rotmat_align()`, `zfs2mat()`, `isstruct()`, `ischar()`.
