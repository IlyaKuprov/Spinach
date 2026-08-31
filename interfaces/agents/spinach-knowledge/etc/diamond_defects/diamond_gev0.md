# etc/diamond_defects/diamond_gev0.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_gev0.m`
- Signature: `[sys,inter]=diamond_gev0(parameters)`
- Total lines: 110

## Purpose

GeV0 spin system for diamond. Syntax: [sys,inter]=diamond_gev0(parameters) Magnetic parameters from Nadolinny et al., Phys. Status Solidi A 213, 2623 (2016), https://doi.org/10.1002/pssa.201600211

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 28-29: Check input count; implemented by `if nargin~=1`.
- Lines 33-34: Check consistency; implemented by `grumble(parameters)`.
- Lines 36-37: Set field-unit conversion constants; implemented by `hz_per_mt=abs(spin('E'))/(2*pi)*1e-3`.
- Lines 39-40: Build the electron tensors; implemented by `electron='E3'`.
- Lines 48-49: Add the germanium isotope if requested; implemented by `germanium=parameters.germanium`.
- Lines 56-57: Build the Spinach structures; implemented by `switch parameters.orientation`.

### Control flow inferred from the code

- Line 29: conditional branch on `nargin~=1`.
- Line 50: conditional branch on `strcmp(germanium,'73Ge')`.
- Line 57: dispatches on `parameters.orientation`; cases `'111'`, `'110'`, `'100'`.
- Line 68: `for` loop over `n=1:numel(nuclei)`.
- Line 76: `for` loop over `n=1:numel(nuclei)`.
- Line 77: conditional branch on `isfield(nuclei{n},'A')&&norm(nuclei{n}.A,2)>0`.

### Key state/data transformations

- Lines 37: computes `hz_per_mt` using `hz_per_mt=abs(spin('E'))/(2*pi)*1e-3`.
- Lines 40: computes `electron` using `electron='E3'`.
- Lines 41-43: computes `frame` using `frame=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 44: computes `gmat` using `gmat=((frame)*diag([2.0027 2.0027 2.0025])*(frame)')`.
- Lines 45: computes `zfs` using `zfs=frame*zfs2mat(80.3*hz_per_mt,0,0,0,0)*frame'`.
- Lines 46: computes `nuclei` using `nuclei={}`.
- Lines 49: computes `germanium` using `germanium=parameters.germanium`.
- Lines 51: computes `nuclei{end+1}` using `nuclei{end+1}=struct('iso','73Ge','A',eye(3)*1.64*hz_per_mt)`.
- Lines 59: computes `C` using `C=rotmat_align([1 1 1],[0 0 1])`.
- Lines 67: computes `sys.isotopes` using `sys.isotopes={electron}`.
- Lines 69: computes `sys.isotopes{n+1}` using `sys.isotopes{n+1}=nuclei{n}.iso`.
- Lines 71: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,numel(sys.isotopes))`.
- Lines 72: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=C*gmat*C'`.
- Lines 73: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes))`.
- Lines 74: computes `[~,~,zfs]` using `[~,~,zfs]=mat2ias(C*zfs*C')`.
- Lines 75: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=zfs`.
- Lines 78: computes `inter.coupling.matrix{1,n+1}` using `inter.coupling.matrix{1,n+1}=C*nuclei{n}.A*C'`.

### Local helper functions

- Line 85: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if(~isstruct(parameters))`.
  - Representative operation: `error('parameters must be a structure.')`.

## Parameters / inputs

- parameters is a structure with the following fields:
- .germanium -'73Ge', 'none', or another germanium isotope
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- GeV0 spin system for diamond. Syntax:
- [sys,inter]=diamond_gev0(parameters)
- Magnetic parameters from Nadolinny et al., Phys. Status Solidi A
- 213, 2623 (2016), https://doi.org/10.1002/pssa.201600211
- parameters is a structure with the following fields:
- .germanium -'73Ge', 'none', or another germanium isotope
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
- Check input count
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `zfs2mat()`, `strcmp()`, `rotmat_align()`, `mat2ias()`, `isfield()`, `isstruct()`, `ischar()`, `ismember()`.
