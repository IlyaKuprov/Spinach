# etc/diamond_defects/diamond_r2.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_r2.m`
- Signature: `[sys,inter]=diamond_r2(parameters)`
- Total lines: 91

## Purpose

R2 self-interstitial spin system for diamond. Syntax: [sys,inter]=diamond_r2(parameters) Magnetic parameters from Hunt et al., Phys. Rev. B 61, 3863 (2000), https://doi.org/10.1103/PhysRevB.61.3863

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check input count; implemented by `if nargin~=1`.
- Lines 33-34: Check consistency; implemented by `grumble(parameters)`.
- Lines 36-37: Build the electron tensors; implemented by `electron='E3'`.
- Lines 44-45: Build the Spinach structures; implemented by `switch parameters.orientation`.

### Control flow inferred from the code

- Line 29: conditional branch on `nargin~=1`.
- Line 45: dispatches on `parameters.orientation`; cases `'111'`, `'110'`, `'100'`.

### Key state/data transformations

- Lines 37: computes `electron` using `electron='E3'`.
- Lines 38-40: computes `frame` using `frame=[0 0 1; 1 0 0; 0 1 0]`.
- Lines 41: computes `gmat` using `gmat=((frame)*diag([2.0019 2.0019 2.0021])*(frame)')`.
- Lines 42: computes `zfs` using `zfs=frame*zfs2mat(parameters.d_sign*4173e6,0,0,0,0)*frame'`.
- Lines 47: computes `C` using `C=rotmat_align([1 1 1],[0 0 1])`.
- Lines 55: computes `sys.isotopes` using `sys.isotopes={electron}`.
- Lines 56: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,numel(sys.isotopes))`.
- Lines 57: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=C*gmat*C'`.
- Lines 58: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes))`.
- Lines 59: computes `[~,~,zfs]` using `[~,~,zfs]=mat2ias(C*zfs*C')`.
- Lines 60: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=zfs`.

### Local helper functions

- Line 65: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if(~isstruct(parameters))`.
  - Representative operation: `error('parameters must be a structure.')`.

## Parameters / inputs

- parameters is a structure with the following fields:
- .d_sign -sign of D
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- R2 self-interstitial spin system for diamond. Syntax:
- [sys,inter]=diamond_r2(parameters)
- Magnetic parameters from Hunt et al., Phys. Rev. B 61,
- 3863 (2000), https://doi.org/10.1103/PhysRevB.61.3863
- parameters is a structure with the following fields:
- .d_sign -sign of D
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
- Check input count
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `zfs2mat()`, `rotmat_align()`, `mat2ias()`, `isstruct()`, `isfield()`, `ischar()`, `any()`, `strcmp()`, `isscalar()`.
