# etc/diamond_defects/diamond_siv0.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_siv0.m`
- Signature: `[sys,inter]=diamond_siv0(parameters)`
- Total lines: 126

## Purpose

SiV0 spin system for diamond. Syntax: [sys,inter]=diamond_siv0(parameters) Magnetic parameters from Edmonds et al., Phys. Rev. B 77, 245205 (2008), https://doi.org/10.1103/PhysRevB.77.245205

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check input count; implemented by `if nargin~=1`.
- Lines 35-36: Check consistency; implemented by `grumble(parameters)`.
- Lines 38-39: Build the electron tensors; implemented by `electron='E3'`.
- Lines 47-48: Add the silicon isotope if requested; implemented by `silicon=parameters.silicon`.
- Lines 55-56: Add reported nearest-neighbour carbons; implemented by `if parameters.n_13c>0`.
- Lines 65-66: Build the Spinach structures; implemented by `switch parameters.orientation`.

### Control flow inferred from the code

- Line 31: conditional branch on `nargin~=1`.
- Line 49: conditional branch on `strcmp(silicon,'29Si')`.
- Line 56: conditional branch on `parameters.n_13c>0`.
- Line 60: `for` loop over `n=1:parameters.n_13c`.
- Line 66: dispatches on `parameters.orientation`; cases `'111'`, `'110'`, `'100'`.
- Line 77: `for` loop over `n=1:numel(nuclei)`.
- Line 85: `for` loop over `n=1:numel(nuclei)`.
- Line 86: conditional branch on `isfield(nuclei{n},'A')&&norm(nuclei{n}.A,2)>0`.

### Key state/data transformations

- Lines 39: computes `electron` using `electron='E3'`.
- Lines 40-42: computes `frame` using `frame=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 43: computes `gmat` using `gmat=((frame)*diag([2.0035 2.0035 2.0042])*(frame)')`.
- Lines 44: computes `zfs` using `zfs=frame*zfs2mat(1000e6,0,0,0,0)*frame'`.
- Lines 45: computes `nuclei` using `nuclei={}`.
- Lines 48: computes `silicon` using `silicon=parameters.silicon`.
- Lines 50: computes `nuclei{end+1}` using `nuclei{end+1}=struct('iso','29Si','A',((frame)*diag([78.9e6 78.9e6 76.3e6])*(frame)'))`.
- Lines 57: computes `Cmat` using `Cmat=((frame)*diag([30.2e6 30.2e6 66.2e6])*(frame)')`.
- Lines 58: computes `nuc_idx` using `nuc_idx=numel(nuclei)`.
- Lines 59: computes `nuclei(nuc_idx+1:nuc_idx+parameters.n_13c)` using `nuclei(nuc_idx+1:nuc_idx+parameters.n_13c)={[]}`.
- Lines 61: computes `nuclei{nuc_idx+n}` using `nuclei{nuc_idx+n}=struct('iso','13C','A',Cmat)`.
- Lines 68: computes `C` using `C=rotmat_align([1 1 1],[0 0 1])`.
- Lines 76: computes `sys.isotopes` using `sys.isotopes={electron}`.
- Lines 78: computes `sys.isotopes{n+1}` using `sys.isotopes{n+1}=nuclei{n}.iso`.
- Lines 80: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,numel(sys.isotopes))`.
- Lines 81: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=C*gmat*C'`.
- Lines 82: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes))`.
- Lines 83: computes `[~,~,zfs]` using `[~,~,zfs]=mat2ias(C*zfs*C')`.

### Local helper functions

- Line 94: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if(~isstruct(parameters))`.
  - Representative operation: `error('parameters must be a structure.')`.

## Parameters / inputs

- parameters is a structure with the following required fields:
- .silicon -'29Si', 'none', or another silicon isotope
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .n_13c -number of reported nearest-neighbour 13C
- hyperfine couplings, between 0 and 6

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- SiV0 spin system for diamond. Syntax:
- [sys,inter]=diamond_siv0(parameters)
- Magnetic parameters from Edmonds et al., Phys. Rev. B 77,
- 245205 (2008), https://doi.org/10.1103/PhysRevB.77.245205
- parameters is a structure with the following required fields:
- .silicon -'29Si', 'none', or another silicon isotope
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .n_13c -number of reported nearest-neighbour 13C
- hyperfine couplings, between 0 and 6
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `zfs2mat()`, `strcmp()`, `nuclei()`, `rotmat_align()`, `mat2ias()`, `isfield()`, `isstruct()`, `ischar()`, `ismember()`, `isscalar()`.
