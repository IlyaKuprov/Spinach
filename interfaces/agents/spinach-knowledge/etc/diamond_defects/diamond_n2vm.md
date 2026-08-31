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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check input count; implemented by `if nargin~=1`.
- Lines 35-36: Check consistency; implemented by `grumble(parameters)`.
- Lines 38-39: Build the reported frames; implemented by `c2rot=anax2dcm([0 0 1],pi)`.
- Lines 44-45: Orthogonalise the carbon hyperfine frame axes; implemented by `xaxis=[-1 -1 0]'`.
- Lines 54-55: Build the electron tensors; implemented by `electron='E'`.
- Lines 59-60: Add nitrogen tensors; implemented by `if strcmp(parameters.nitrogen,'15N')`.
- Lines 76-77: Add reported nearest-neighbour carbons; implemented by `if parameters.include_13c`.
- Lines 82-83: Build the Spinach structures; implemented by `switch parameters.orientation`.

### Control flow inferred from the code

- Line 31: conditional branch on `nargin~=1`.
- Line 60: conditional branch on `strcmp(parameters.nitrogen,'15N')`.
- Line 77: conditional branch on `parameters.include_13c`.
- Line 83: dispatches on `parameters.orientation`; cases `'111'`, `'110'`, `'100'`.
- Line 94: `for` loop over `n=1:numel(nuclei)`.
- Line 100: `for` loop over `n=1:numel(nuclei)`.
- Line 101: conditional branch on `isfield(nuclei{n},'A')&&norm(nuclei{n}.A,2)>0`.
- Line 104: conditional branch on `isfield(nuclei{n},'Q')&&~isempty(nuclei{n}.Q)`.

### Key state/data transformations

- Lines 39: computes `c2rot` using `c2rot=anax2dcm([0 0 1],pi)`.
- Lines 40: computes `gframe` using `gframe=diamond_frame_xyz([1 1 0],[0 0 1],[1 -1 0])`.
- Lines 41-42: computes `nframe` using `nframe=anax2dcm([1 -1 0],-3.5*pi/180)* diamond_frame_xyz([1 1 -2],[1 1 1],[1 -1 0])`.
- Lines 45: computes `xaxis` using `xaxis=[-1 -1 0]'`.
- Lines 47: computes `zaxis` using `zaxis=[-1 1 1]'`.
- Lines 50: computes `yaxis` using `yaxis=cross(zaxis,xaxis)`.
- Lines 51-52: computes `cframe` using `cframe=anax2dcm([-1 -1 0],2.0*pi/180)* diamond_frame_xyz(xaxis,yaxis,zaxis)`.
- Lines 55: computes `electron` using `electron='E'`.
- Lines 56: computes `gmat` using `gmat=((gframe)*diag([2.00345 2.00274 2.00271])*(gframe)')`.
- Lines 57: computes `nuclei` using `nuclei={}`.
- Lines 61: computes `nuclei{end+1}` using `nuclei{end+1}=struct('iso','15N','A',((nframe)*diag([3.47e6 4.51e6 4.09e6])*(nframe)'))`.
- Lines 64: computes `scale` using `scale=spin('14N')/spin('15N')`.
- Lines 65-67: computes `Qframe` using `Qframe=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 85: computes `C` using `C=rotmat_align([1 1 1],[0 0 1])`.
- Lines 93: computes `sys.isotopes` using `sys.isotopes={electron}`.
- Lines 95: computes `sys.isotopes{n+1}` using `sys.isotopes{n+1}=nuclei{n}.iso`.
- Lines 97: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,numel(sys.isotopes))`.
- Lines 98: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=C*gmat*C'`.

### Local helper functions

- Line 113: `diamond_frame_xyz()` — `function frame=diamond_frame_xyz(xaxis,yaxis,zaxis)`. Consistency enforcement
  - Representative operation: `frame=[xaxis(:) yaxis(:) zaxis(:)]`.
  - Representative operation: `[frame,~]=qr(frame,0)`.
- Line 122: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if(~isstruct(parameters))`.
  - Representative operation: `error('parameters must be a structure.')`.

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
