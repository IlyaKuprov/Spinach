# etc/diamond_defects/diamond_n_inter.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_n_inter.m`
- Signature: `[sys,inter]=diamond_n_inter(parameters)`
- Total lines: 131

## Purpose

Nitrogen interstitial spin system for diamond. Syntax: [sys,inter]=diamond_n_inter(parameters) Magnetic parameters from Felton et al., J. Phys. Condens. Matter 21, 364212 (2009), https://doi.org/10.1088/0953-8984/21/36/364212

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `diamond_frame_xyz()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 29-30: Check input count; implemented by `if nargin~=1`.
- Lines 34-35: Check consistency; implemented by `grumble(parameters)`.
- Lines 37-38: Build the common frame; implemented by `electron='E'`.
- Lines 43-44: Select the interstitial centre; implemented by `centre=lower(parameters.centre)`.
- Lines 59-60: Add the nitrogen hyperfine tensor; implemented by `switch parameters.nitrogen`.
- Lines 70-71: Build the Spinach structures; implemented by `switch parameters.orientation`.

### Control flow inferred from the code

- Line 30: conditional branch on `nargin~=1`.
- Line 45: dispatches on `centre`; cases `'war9'`, `'war10'`.
- Line 60: dispatches on `parameters.nitrogen`; cases `'15N'`, `'14N'`.
- Line 71: dispatches on `parameters.orientation`; cases `'111'`, `'110'`, `'100'`.

### Key state/data transformations

- Lines 38: computes `electron` using `electron='E'`.
- Lines 39-41: computes `gframe` using `gframe=diamond_frame_xyz([sind(90)*cosd(45);sind(90)*sind(45);cosd(90)], [sind(180)*cosd(45);sind(180)*sind(45);cosd(180)], [sind(90)*cosd(315);sind(90)*sind(315);cosd(9…`.
- Lines 44: computes `centre` using `centre=lower(parameters.centre)`.
- Lines 47: computes `gmat` using `gmat=((gframe)*diag([2.00343 2.00272 2.00268])*(gframe)')`.
- Lines 48: computes `Amat` using `Amat=((gframe)*diag([8.30e6 7.85e6 8.17e6])*(gframe)')`.
- Lines 50-52: computes `aframe` using `aframe=diamond_frame_xyz([sind(44.8)*cosd(45.0);sind(44.8)*sind(45.0);cosd(44.8)], [sind(134.8)*cosd(45.0);sind(134.8)*sind(45.0);cosd(134.8)], [sind(90)*cosd(315);sind(…`.
- Lines 62: computes `nucleus` using `nucleus=struct('iso','15N','A',Amat)`.
- Lines 64: computes `scale` using `scale=spin('14N')/spin('15N')`.
- Lines 73: computes `C` using `C=rotmat_align([1 1 1],[0 0 1])`.
- Lines 81: computes `sys.isotopes` using `sys.isotopes={electron,nucleus.iso}`.
- Lines 82: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,numel(sys.isotopes))`.
- Lines 83: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=C*gmat*C'`.
- Lines 84: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes))`.
- Lines 85: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=C*nucleus.A*C'`.

### Local helper functions

- Line 90: `diamond_frame_xyz()` — `function frame=diamond_frame_xyz(xaxis,yaxis,zaxis)`. Consistency enforcement
  - Representative operation: `frame=[xaxis(:) yaxis(:) zaxis(:)]`.
  - Representative operation: `[frame,~]=qr(frame,0)`.
- Line 99: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if(~isstruct(parameters))`.
  - Representative operation: `error('parameters must be a structure.')`.

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'war9' or 'war10'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nitrogen -'14N' or '15N'

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Nitrogen interstitial spin system for diamond. Syntax:
- [sys,inter]=diamond_n_inter(parameters)
- Magnetic parameters from Felton et al., J. Phys. Condens. Matter
- 21, 364212 (2009), https://doi.org/10.1088/0953-8984/21/36/364212
- parameters is a structure with the following fields:
- .centre -'war9' or 'war10'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nitrogen -'14N' or '15N'
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
- Check input count

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `diamond_frame_xyz()`, `sind()`, `cosd()`, `lower()`, `spin()`, `rotmat_align()`, `xaxis()`, `yaxis()`, `zaxis()`, `frame()`, `isstruct()`, `isfield()`, `ischar()`, `any()`, `strcmp()`.
