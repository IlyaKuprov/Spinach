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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 37-38: Check input count; implemented by `if nargin~=1`.
- Lines 42-43: Check consistency; implemented by `grumble(parameters)`.
- Lines 45-46: Select the vacancy centre; implemented by `centre=lower(parameters.centre)`.
- Lines 48-49: Set R4/W6 principal-axis frame; implemented by `r4_x=[1;-1;0]/sqrt(2)`.
- Lines 56-57: Set W29 principal-axis frame; implemented by `w29_x=[0;-1;1]/sqrt(2)`.
- Lines 64-67: Set vacancy-chain principal-axis frame; implemented by `frame_110=[-1/sqrt(2) 0 1/sqrt(2); 1/sqrt(2) 0 1/sqrt(2); 0 1 0]`.
- Lines 101-102: Build the Spinach structures; implemented by `switch parameters.orientation`.

### Control flow inferred from the code

- Line 38: conditional branch on `nargin~=1`.
- Line 68: dispatches on `centre`; cases `{'r4_w6','w6','r4'}`, `'w29'`, `'r5'`, `'o1'`, `'r6'`, `'r10'`, `'r11'`.
- Line 102: dispatches on `parameters.orientation`; cases `'111'`, `'110'`, `'100'`.

### Key state/data transformations

- Lines 46: computes `centre` using `centre=lower(parameters.centre)`.
- Lines 49: computes `r4_x` using `r4_x=[1;-1;0]/sqrt(2)`.
- Lines 50: computes `r4_z` using `r4_z=[sind(54.2)*cosd(45);sind(54.2)*sind(45);cosd(54.2)]`.
- Lines 52: computes `r4_y` using `r4_y=cross(r4_z,r4_x)`.
- Lines 54: computes `frame_r4` using `frame_r4=[r4_x r4_y r4_z]`.
- Lines 57: computes `w29_x` using `w29_x=[0;-1;1]/sqrt(2)`.
- Lines 58: computes `w29_z` using `w29_z=[0.619;-0.556;-0.556]`.
- Lines 60: computes `w29_y` using `w29_y=cross(w29_z,w29_x)`.
- Lines 62: computes `frame_w29` using `frame_w29=[w29_x w29_y w29_z]`.
- Lines 65-67: computes `frame_110` using `frame_110=[-1/sqrt(2) 0 1/sqrt(2); 1/sqrt(2) 0 1/sqrt(2); 0 1 0]`.
- Lines 70: computes `electron` using `electron='E3'`.
- Lines 71: computes `gmat` using `gmat=((frame_r4)*diag([2.0022 2.0026 2.0013])*(frame_r4)')`.
- Lines 72: computes `zfs` using `zfs=((frame_r4)*diag([105 197 -303]*1e6)*(frame_r4)')`.
- Lines 104: computes `C` using `C=rotmat_align([1 1 1],[0 0 1])`.
- Lines 112: computes `sys.isotopes` using `sys.isotopes={electron}`.
- Lines 113: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,numel(sys.isotopes))`.
- Lines 114: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=C*gmat*C'`.
- Lines 115: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes))`.

### Local helper functions

- Line 122: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if(~isstruct(parameters))`.
  - Representative operation: `error('parameters must be a structure.')`.

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
