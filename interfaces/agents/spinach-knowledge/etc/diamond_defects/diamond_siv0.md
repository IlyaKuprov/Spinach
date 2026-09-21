# etc/diamond_defects/diamond_siv0.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_siv0.m`
- Signature: `[sys,inter]=diamond_siv0(parameters)`
- Total lines: 126

## Purpose

SiV0 spin system for diamond. Syntax: [sys,inter]=diamond_siv0(parameters) Magnetic parameters from Edmonds et al., Phys. Rev. B 77, 245205 (2008), https://doi.org/10.1103/PhysRevB.77.245205

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters is a structure with the following required fields:
- .silicon -'29Si', 'none', or another silicon isotope
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .n_13c -number of reported nearest-neighbour 13C
- hyperfine couplings, between 0 and 6

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- SiV0 spin system for diamond. Syntax:
- [sys,inter]=diamond_siv0(parameters)
- Magnetic parameters from Edmonds et al., Phys. Rev. B 77,
- 245205 (2008), https://doi.org/10.1103/PhysRevB.77.245205
- parameters is a structure with the following required fields:
- .silicon -'29Si', 'none', or another silicon isotope
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .n_13c -number of reported nearest-neighbour 13C
- hyperfine couplings, between 0 and 6
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `zfs2mat()`, `strcmp()`, `nuclei()`, `rotmat_align()`, `mat2ias()`, `isfield()`, `isstruct()`, `ischar()`, `ismember()`, `isscalar()`.
