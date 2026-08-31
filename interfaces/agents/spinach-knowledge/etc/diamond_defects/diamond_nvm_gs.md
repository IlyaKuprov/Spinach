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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(parameters)`.
- Lines 35-36: Set default parameters; implemented by `if ~isfield(parameters,'orientation')`.
- Lines 43-44: Define spin system isotopes; implemented by `sys.isotopes={'E3',parameters.nitrogen}`.
- Lines 48-51: Set the trigonal principal-axis frame; implemented by `frame=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 53-54: Rotation matrix for orientation; implemented by `switch parameters.orientation`.
- Lines 70-71: Complain and bomb out; implemented by `error('unknown orientation specification')`.
- Lines 75-76: Electron g-tensor; implemented by `inter.zeeman.matrix{1}=C*frame*diag([2.0031 2.0031 2.0029])*frame'*C'`.
- Lines 78-79: Electron ZFS tensor; implemented by `inter.coupling.matrix{1,1}=C*frame*zfs2mat(2872e6,0,0,0,0)*frame'*C'`.
- Lines 81-82: HFC and NQI tensor; implemented by `switch parameters.nitrogen`.
- Lines 86-87: Hyperfine and quadrupolar coupling; implemented by `inter.coupling.matrix{1,2}=C*frame*diag([-2.70e6 -2.70e6 -2.14e6])*frame'*C'`.
- Lines 92-93: Only hyperfine coupling, opposite sign; implemented by `inter.coupling.matrix{1,2}=C*frame*diag([3.65e6 3.65e6 3.03e6])*frame'*C'`.
- Lines 97-98: Complain and bomb out; implemented by `error('wrong nitrogen isotope.')`.

### Control flow inferred from the code

- Line 36: conditional branch on `~isfield(parameters,'orientation')`.
- Line 39: conditional branch on `~isfield(parameters,'nitrogen')`.
- Line 54: dispatches on `parameters.orientation`; cases `'100'`, `'110'`, `'111'`.
- Line 82: dispatches on `parameters.nitrogen`; cases `'14N'`, `'15N'`.

### Key state/data transformations

- Lines 37: computes `parameters.orientation` using `parameters.orientation='111'`.
- Lines 40: computes `parameters.nitrogen` using `parameters.nitrogen='14N'`.
- Lines 44: computes `sys.isotopes` using `sys.isotopes={'E3',parameters.nitrogen}`.
- Lines 45: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,numel(sys.isotopes))`.
- Lines 46: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes))`.
- Lines 49-51: computes `frame` using `frame=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 58: computes `C` using `C=rotmat_align([1 0 0],[0 0 1])`.
- Lines 76: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=C*frame*diag([2.0031 2.0031 2.0029])*frame'*C'`.
- Lines 79: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=C*frame*zfs2mat(2872e6,0,0,0,0)*frame'*C'`.
- Lines 87: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=C*frame*diag([-2.70e6 -2.70e6 -2.14e6])*frame'*C'`.
- Lines 88: computes `inter.coupling.matrix{2,2}` using `inter.coupling.matrix{2,2}=C*frame*zfs2mat(-5.01e6,0,0,0,0)*frame'*C'`.

### Local helper functions

- Line 105: `grumble()` — `function grumble(parameters)`. Once, when sheltering under a tree during a storm near Lichfield,
  - Representative operation: `if(~isstruct(parameters))`.
  - Representative operation: `error('parameters must be a structure.')`.

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
