# etc/diamond_defects/diamond_p.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_p.m`
- Signature: `[sys,inter]=diamond_p(parameters)`
- Total lines: 170

## Purpose

Phosphorus-related defect spin system for diamond. Syntax: [sys,inter]=diamond_p(parameters) Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check input count; implemented by `if nargin~=1`.
- Lines 36-37: Set carbon coupling default; implemented by `if isstruct(parameters)&&~isfield(parameters,'include_13c')`.
- Lines 41-42: Check consistency; implemented by `grumble(parameters)`.
- Lines 44-45: Set field-unit conversion constants; implemented by `hz_per_mt=abs(spin('E'))/(2*pi)*1e-3`.
- Lines 47-48: Select the phosphorus centre; implemented by `centre=lower(parameters.centre)`.
- Lines 52-55: Set the trigonal principal-axis frame; implemented by `frame_111=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 57-58: Get phosphorus centre table data; implemented by `switch centre`.
- Lines 107-108: Build the Spinach structures; implemented by `switch parameters.orientation`.

### Control flow inferred from the code

- Line 32: conditional branch on `nargin~=1`.
- Line 37: conditional branch on `isstruct(parameters)&&~isfield(parameters,'include_13c')`.
- Line 58: dispatches on `centre`; cases `'ma1'`, `'np1'`, `'np2'`, `'np3'`, `'np4'`, `'np5'`, `'np6'`, `'np8'`, `'np9'`.
- Line 63: conditional branch on `parameters.include_13c`.
- Line 108: dispatches on `parameters.orientation`; cases `'111'`, `'110'`, `'100'`.
- Line 119: `for` loop over `n=1:numel(nuclei)`.
- Line 125: `for` loop over `n=1:numel(nuclei)`.
- Line 126: conditional branch on `isfield(nuclei{n},'A')&&norm(nuclei{n}.A,2)>0`.

### Key state/data transformations

- Lines 38: computes `parameters.include_13c` using `parameters.include_13c=false`.
- Lines 45: computes `hz_per_mt` using `hz_per_mt=abs(spin('E'))/(2*pi)*1e-3`.
- Lines 48: computes `centre` using `centre=lower(parameters.centre)`.
- Lines 49: computes `electron` using `electron='E'`.
- Lines 50: computes `nuclei` using `nuclei={}; frame=eye(3)`.
- Lines 53-55: computes `frame_111` using `frame_111=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 60: computes `gmat` using `gmat=eye(3)*2.0025`.
- Lines 61-62: computes `nuclei{end+1}` using `nuclei{end+1}=struct('iso','31P','A', ((frame_111)*diag([1.96 1.96 2.32]*hz_per_mt)*(frame_111)'))`.
- Lines 110: computes `C` using `C=rotmat_align([1 1 1],[0 0 1])`.
- Lines 118: computes `sys.isotopes` using `sys.isotopes={electron}`.
- Lines 120: computes `sys.isotopes{n+1}` using `sys.isotopes{n+1}=nuclei{n}.iso`.
- Lines 122: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,numel(sys.isotopes))`.
- Lines 123: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=C*gmat*C'`.
- Lines 124: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes))`.
- Lines 127: computes `inter.coupling.matrix{1,n+1}` using `inter.coupling.matrix{1,n+1}=C*nuclei{n}.A*C'`.

### Local helper functions

- Line 134: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if(~isstruct(parameters))`.
  - Representative operation: `error('parameters must be a structure.')`.

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'ma1', 'np1', 'np2', 'np3', 'np4', 'np5',
- 'np6', 'np8', or 'np9'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .include_13c -include the reported 13C hyperfine coupling;
- applies only to MA1 and defaults to false

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Phosphorus-related defect spin system for diamond. Syntax:
- [sys,inter]=diamond_p(parameters)
- Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),
- parameters is a structure with the following fields:
- .centre -'ma1', 'np1', 'np2', 'np3', 'np4', 'np5',
- 'np6', 'np8', or 'np9'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .include_13c -include the reported 13C hyperfine coupling;
- applies only to MA1 and defaults to false
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `isstruct()`, `isfield()`, `grumble()`, `spin()`, `lower()`, `rotmat_align()`, `ischar()`, `ismember()`, `islogical()`, `isscalar()`, `strcmpi()`.
