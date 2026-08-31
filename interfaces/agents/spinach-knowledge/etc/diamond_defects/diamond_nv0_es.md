# etc/diamond_defects/diamond_nv0_es.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_nv0_es.m`
- Signature: `[sys,inter]=diamond_nv0_es(parameters)`
- Total lines: 108

## Purpose

NV0 excited-state spin system for diamond. Syntax: [sys,inter]=diamond_nv0_es(parameters) Magnetic parameters from Felton et al., Phys. Rev. B 77, 081201 (2008), https://doi.org/10.1103/PhysRevB.77.081201

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check input count; implemented by `if nargin~=1`.
- Lines 35-36: Check consistency; implemented by `grumble(parameters)`.
- Lines 38-41: Build the electron tensors; implemented by `frame=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 46-47: Add the nitrogen hyperfine tensor; implemented by `Amat=((frame)*diag([-23.8e6 -23.8e6 -35.7e6])*(frame)')`.
- Lines 58-59: Build the Spinach structures; implemented by `switch parameters.orientation`.

### Control flow inferred from the code

- Line 31: conditional branch on `nargin~=1`.
- Line 48: dispatches on `parameters.nitrogen`; cases `'15N'`, `'14N'`.
- Line 59: dispatches on `parameters.orientation`; cases `'111'`, `'110'`, `'100'`.

### Key state/data transformations

- Lines 39-41: computes `frame` using `frame=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 42: computes `electron` using `electron='E4'`.
- Lines 43: computes `gmat` using `gmat=((frame)*diag([2.0035 2.0035 2.0029])*(frame)')`.
- Lines 44: computes `zfs` using `zfs=frame*zfs2mat(1685e6,0,0,0,0)*frame'`.
- Lines 47: computes `Amat` using `Amat=((frame)*diag([-23.8e6 -23.8e6 -35.7e6])*(frame)')`.
- Lines 50: computes `nucleus` using `nucleus=struct('iso','15N','A',Amat)`.
- Lines 52: computes `scale` using `scale=spin('14N')/spin('15N')`.
- Lines 61: computes `C` using `C=rotmat_align([1 1 1],[0 0 1])`.
- Lines 69: computes `sys.isotopes` using `sys.isotopes={electron,nucleus.iso}`.
- Lines 70: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,numel(sys.isotopes))`.
- Lines 71: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=C*gmat*C'`.
- Lines 72: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes))`.
- Lines 73: computes `[~,~,zfs]` using `[~,~,zfs]=mat2ias(C*zfs*C')`.
- Lines 74: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=zfs`.
- Lines 75: computes `inter.coupling.matrix{1,2}` using `inter.coupling.matrix{1,2}=C*nucleus.A*C'`.

### Local helper functions

- Line 80: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if(~isstruct(parameters))`.
  - Representative operation: `error('parameters must be a structure.')`.

## Parameters / inputs

- parameters is a structure with the following required fields:
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nitrogen -'14N' or '15N'; 14N hyperfine couplings
- are scaled from 15N, with no NQI included
- because none is reported for this state

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- NV0 excited-state spin system for diamond. Syntax:
- [sys,inter]=diamond_nv0_es(parameters)
- Magnetic parameters from Felton et al., Phys. Rev. B 77,
- 081201 (2008), https://doi.org/10.1103/PhysRevB.77.081201
- parameters is a structure with the following required fields:
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .nitrogen -'14N' or '15N'; 14N hyperfine couplings
- are scaled from 15N, with no NQI included
- because none is reported for this state
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `zfs2mat()`, `spin()`, `rotmat_align()`, `mat2ias()`, `isstruct()`, `isfield()`, `ischar()`, `ismember()`, `any()`, `strcmp()`.
