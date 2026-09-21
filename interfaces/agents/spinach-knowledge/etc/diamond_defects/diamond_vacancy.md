# etc/diamond_defects/diamond_vacancy.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_vacancy.m`
- Signature: `[sys,inter]=diamond_vacancy(parameters)`
- Total lines: 144

## Purpose

Vacancy-family defect spin systems for diamond. Syntax: [sys,inter]=diamond_vacancy(parameters) R4/W6 parameters from Twitchen et al., Phys. Rev. B 59, 12900 (1999), https://doi.org/10.1103/PhysRevB.59.12900; W29 parameters from Kirui et al., Diam. Relat. Mater. 8, 1569 (1999), https://doi.org/10.1016/S0925-9635(99)00037-0; R5/O1/R6/R10/R11 parameters from Iakoubovskii and Stesmans, Phys. Rev. B 66, 045406 (2002), cr

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'r4_w6', 'w6', 'r4', 'w29', 'r5', 'o1',
- 'r6', 'r10', or 'r11'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Vacancy-family defect spin systems for diamond. Syntax:
- [sys,inter]=diamond_vacancy(parameters)
- R4/W6 parameters from Twitchen et al., Phys. Rev. B 59, 12900
- (1999), https://doi.org/10.1103/PhysRevB.59.12900; W29
- parameters from Kirui et al., Diam. Relat. Mater. 8, 1569
- (1999), https://doi.org/10.1016/S0925-9635(99)00037-0;
- R5/O1/R6/R10/R11 parameters from Iakoubovskii and Stesmans,
- Phys. Rev. B 66, 045406 (2002),
- cross-checked against Ball, PhD thesis, OIST Graduate University
- (2021).
- parameters is a structure with the following fields:
- .centre -'r4_w6', 'w6', 'r4', 'w29', 'r5', 'o1',

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `lower()`, `sind()`, `cosd()`, `cross()`, `rotmat_align()`, `mat2ias()`, `isstruct()`, `isfield()`, `ischar()`.
