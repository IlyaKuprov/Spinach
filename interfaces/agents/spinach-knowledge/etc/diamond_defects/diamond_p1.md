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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(parameters)`.
- Lines 35-36: Set default parameters; implemented by `if ~isfield(parameters,'orientation')`.
- Lines 43-44: Define spin system isotopes; implemented by `sys.isotopes={'E',parameters.nitrogen}`.
- Lines 48-51: Set the trigonal principal-axis frame; implemented by `frame=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 53-54: Rotation matrix for orientation; implemented by `switch parameters.orientation`.
- Lines 70-71: Complain and bomb out; implemented by `error('unknown oritentation specification')`.
- Lines 75-76: Electron g-tensor; implemented by `inter.zeeman.matrix{1}=C*frame*diag([2.00220 2.00220 2.00218])*frame'*C'`.
- Lines 78-79: HFC and NQI tensor; implemented by `switch parameters.nitrogen`.
- Lines 83-84: Hyperfine and quadrupolar coupling; implemented by `inter.coupling.matrix{1,2}=C*frame*diag([81.3e6 81.3e6 114.0e6])*frame'*C'`.
- Lines 89-90: Only hyperfine coupling, opposite sign; implemented by `inter.coupling.matrix{1,2}=C*frame*diag([-114.0e6 -114.0e6 -159.9e6])*frame'*C'`.
- Lines 94-95: Complain and bomb out; implemented by `error('wrong nitrogen isotope.')`.

### Control flow inferred from the code

- Line 36: conditional branch on `~isfield(parameters,'orientation')`.
- Line 39: conditional branch on `~isfield(parameters,'nitrogen')`.
- Line 54: dispatches on `parameters.orientation`; cases `'100'`, `'110'`, `'111'`.
- Line 79: dispatches on `parameters.nitrogen`; cases `'14N'`, `'15N'`.

### Key state/data transformations

- Lines 37: computes `parameters.orientation` using `parameters.orientation='111'`.
- Lines 40: computes `parameters.nitrogen` using `parameters.nitrogen='14N'`.
- Lines 44: computes `sys.isotopes` using `sys.isotopes={'E',parameters.nitrogen}`.
- Lines 45: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,numel(sys.isotopes))`.
- Lines 46: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes))`.
- Lines 49-51: computes `frame` using `frame=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 58: computes `C` using `C=rotmat_align([1 0 0],[0 0 1])`.
- Lines 76: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=C*frame*diag([2.00220 2.00220 2.00218])*frame'*C'`.
- Lines 84: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=C*frame*diag([81.3e6 81.3e6 114.0e6])*frame'*C'`.
- Lines 85: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=C*frame*zfs2mat(-3.97e6,0,0,0,0)*frame'*C'`.

### Local helper functions

- Line 102: `grumble()` — `function grumble(parameters)`. One observes the survivors and learns from them.
  - Representative operation: `if(~isstruct(parameters))`.
  - Representative operation: `error('parameters must be a structure.')`.

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
