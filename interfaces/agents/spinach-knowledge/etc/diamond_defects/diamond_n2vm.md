# etc/diamond_defects/diamond_n2vm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_n2vm.m`
- Signature: `[sys,inter]=diamond_n2vm(parameters)`
- Total lines: 155

## Purpose

N2V-spin system for diamond. Syntax: [sys,inter]=diamond_n2vm(parameters) Magnetic parameters from Green et al., Phys. Rev. B 92, 165204 (2015), https://doi.org/10.1103/PhysRevB.92.165204

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `diamond_frame_xyz()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters is a structure with the following required fields:
- .nitrogen -'14N' or '15N'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .include_13c -include reported 13C hyperfine couplings,
- true or false

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- N2V-spin system for diamond. Syntax:
- [sys,inter]=diamond_n2vm(parameters)
- Magnetic parameters from Green et al., Phys. Rev. B 92,
- 165204 (2015), https://doi.org/10.1103/PhysRevB.92.165204
- parameters is a structure with the following required fields:
- .nitrogen -'14N' or '15N'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .include_13c -include reported 13C hyperfine couplings,
- true or false
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `anax2dcm()`, `diamond_frame_xyz()`, `dot()`, `cross()`, `strcmp()`, `spin()`, `zfs2mat()`, `rotmat_align()`, `isfield()`, `mat2ias()`, `xaxis()`, `yaxis()`, `zaxis()`, `frame()`, `isstruct()`.
