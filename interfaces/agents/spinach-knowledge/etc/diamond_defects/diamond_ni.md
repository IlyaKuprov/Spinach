# etc/diamond_defects/diamond_ni.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_ni.m`
- Signature: `[sys,inter]=diamond_ni(parameters)`
- Total lines: 232

## Purpose

Nickel-related defect spin system for diamond. Syntax: [sys,inter]=diamond_ni(parameters) W8 magnetic parameters from Ludwig and Woodbury, Phys. Rev. B 41, 3905 (1990), https://doi.org/10.1103/PhysRevB.41.3905 The W8 quartet entry assumes zero ZFS: no ZFS parameters were reported in the cited data, and off-central transitions are treated as unresolved rather than explicitly modelled. Other nickel-centre table values 

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `diamond_frame_xyz()`, `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 39-40: Check input count; implemented by `if nargin~=1`.
- Lines 44-45: Check consistency; implemented by `grumble(parameters)`.
- Lines 47-48: Set irrelevant carbon count to zero; implemented by `if ~strcmpi(parameters.centre,'w8')`.
- Lines 52-53: Set field-unit conversion constants; implemented by `hz_per_mt=abs(spin('E'))/(2*pi)*1e-3`.
- Lines 56-57: Select the nickel centre; implemented by `centre=lower(parameters.centre)`.
- Lines 59-62: Set the trigonal principal-axis frame; implemented by `frame_111=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 85-86: Get tabulated NE centre parameters; implemented by `switch centre`.
- Lines 104-105: Build the Nadolinny table frame; implemented by `xaxis=[1;-1;0]/sqrt(2)`.
- Lines 123-124: Get tabulated AB centre parameters; implemented by `switch centre`.
- Lines 149-150: Build the Spinach structures; implemented by `switch parameters.orientation`.

### Control flow inferred from the code

- Line 40: conditional branch on `nargin~=1`.
- Line 48: conditional branch on `~strcmpi(parameters.centre,'w8')`.
- Line 64: dispatches on `centre`; cases `'w8'`, `{'ne1','ne2','ne3','ne5','ne8'}`, `'ne1'`, `'ne2'`, `'ne3'`, `'ne5'`, `'ne8'`, `'ne4'`, `{'ab1','ab2','ab3','ab4'}`, `'ab1'`, ….
- Line 69: conditional branch on `strcmp(nickel,'61Ni')`.
- Line 74: conditional branch on `parameters.n_13c>0`.
- Line 78: `for` loop over `n=1:parameters.n_13c`.
- Line 86: dispatches on `centre`; cases `'ne1'`, `'ne2'`, `'ne3'`, `'ne5'`, `'ne8'`, `'ne4'`, `{'ab1','ab2','ab3','ab4'}`, `'ab1'`, `'ab2'`, `'ab3'`, ….
- Line 113: `for` loop over `n=1:size(avalues,1)`.
- Line 124: dispatches on `centre`; cases `'ab1'`, `'ab2'`, `'ab3'`, `'ab4'`, `'ab5'`, `{'nol1','nirim5'}`.
- Line 150: dispatches on `parameters.orientation`; cases `'111'`, `'110'`, `'100'`.
- Line 161: `for` loop over `n=1:numel(nuclei)`.
- Line 167: conditional branch on `~isempty(zfs)`.
- Line 171: `for` loop over `n=1:numel(nuclei)`.
- Line 172: conditional branch on `isfield(nuclei{n},'A')&&norm(nuclei{n}.A,2)>0`.

### Key state/data transformations

- Lines 49: computes `parameters.n_13c` using `parameters.n_13c=0`.
- Lines 53: computes `hz_per_mt` using `hz_per_mt=abs(spin('E'))/(2*pi)*1e-3`.
- Lines 54: computes `hz_per_t` using `hz_per_t=abs(spin('E'))/(2*pi)`.
- Lines 57: computes `centre` using `centre=lower(parameters.centre)`.
- Lines 60-62: computes `frame_111` using `frame_111=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 63: computes `nuclei` using `nuclei={}; zfs=[]`.
- Lines 66: computes `electron` using `electron='E4'`.
- Lines 67: computes `gmat` using `gmat=eye(3)*2.032`.
- Lines 68: computes `nickel` using `nickel=parameters.nickel`.
- Lines 70: computes `nuclei{end+1}` using `nuclei{end+1}=struct('iso','61Ni','A',eye(3)*0.65*hz_per_mt)`.
- Lines 75: computes `Cmat` using `Cmat=((frame_111)*diag([0.340 0.340 1.339]*hz_per_mt)*(frame_111)')`.
- Lines 76: computes `nuc_idx` using `nuc_idx=numel(nuclei)`.
- Lines 77: computes `nuclei(nuc_idx+1:nuc_idx+parameters.n_13c)` using `nuclei(nuc_idx+1:nuc_idx+parameters.n_13c)={[]}`.
- Lines 79: computes `nuclei{nuc_idx+n}` using `nuclei{nuc_idx+n}=struct('iso','13C','A',Cmat)`.
- Lines 88: computes `gvals` using `gvals=[2.1282 2.0070 2.0908]; alpha=14`.
- Lines 89: computes `avalues` using `avalues=[2.09 1.43 1.45;2.09 1.43 1.45]`.
- Lines 105: computes `xaxis` using `xaxis=[1;-1;0]/sqrt(2)`.
- Lines 106: computes `ybase` using `ybase=[1;1;0]/sqrt(2)`.

### Local helper functions

- Line 180: `diamond_frame_xyz()` — `function frame=diamond_frame_xyz(xaxis,yaxis,zaxis)`. Consistency enforcement
  - Representative operation: `frame=[xaxis(:) yaxis(:) zaxis(:)]`.
  - Representative operation: `[frame,~]=qr(frame,0)`.
- Line 189: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if ~isstruct(parameters)`.
  - Representative operation: `error('parameters must be a structure.')`.

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'w8', 'ne1', 'ne2', 'ne3', 'ne4', 'ne5',
- 'ne8', 'ab1', 'ab2', 'ab3', 'ab4', 'ab5',
- 'nol1', or 'nirim5'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nickel -'61Ni', 'none', or another nickel isotope;
- required when .centre is 'w8'
- .n_13c -number of reported 13C hyperfine couplings
- to include, from 0 to 4; applies only to W8

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Nickel-related defect spin system for diamond. Syntax:
- [sys,inter]=diamond_ni(parameters)
- W8 magnetic parameters from Ludwig and Woodbury, Phys. Rev. B 41,
- 3905 (1990), https://doi.org/10.1103/PhysRevB.41.3905
- The W8 quartet entry assumes zero ZFS: no ZFS parameters were
- reported in the cited data, and off-central transitions are treated
- as unresolved rather than explicitly modelled.
- Other nickel-centre table values from Nadolinny et al., Crystals
- 7, 237 (2017), https://doi.org/10.3390/cryst7080237
- parameters is a structure with the following fields:
- .centre -'w8', 'ne1', 'ne2', 'ne3', 'ne4', 'ne5',
- 'ne8', 'ab1', 'ab2', 'ab3', 'ab4', 'ab5',

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmpi()`, `spin()`, `lower()`, `strcmp()`, `nuclei()`, `cosd()`, `sind()`, `cross()`, `diamond_frame_xyz()`, `avalues()`, `zfs2mat()`, `rotmat_align()`, `mat2ias()`, `isfield()`, `xaxis()`.
