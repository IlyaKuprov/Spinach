# etc/diamond_defects/diamond_co.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_co.m`
- Signature: `[sys,inter]=diamond_co(parameters)`
- Total lines: 125

## Purpose

Cobalt-related defect spin system for diamond. Syntax: [sys,inter]=diamond_co(parameters) Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),

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
- Lines 39-40: Select the cobalt centre; implemented by `centre=lower(parameters.centre)`.
- Lines 43-44: Get cobalt centre table data; implemented by `switch centre`.
- Lines 57-58: Assemble the cobalt principal-axis frame; implemented by `yaxis=[0;1;1]/sqrt(2)`.
- Lines 65-66: Orthogonalise the frame; implemented by `[frame,~]=qr(frame,0)`.
- Lines 68-69: Enforce right-handed axes; implemented by `if det(frame)<0`.
- Lines 73-74: Build electron and cobalt tensors; implemented by `gmat=((frame)*diag(gvals)*(frame)')`.
- Lines 77-78: Build the Spinach structures; implemented by `switch parameters.orientation`.

### Control flow inferred from the code

- Line 29: conditional branch on `nargin~=1`.
- Line 44: dispatches on `centre`; cases `'o4'`, `'nlo2'`.
- Line 69: conditional branch on `det(frame)<0`.
- Line 78: dispatches on `parameters.orientation`; cases `'111'`, `'110'`, `'100'`.

### Key state/data transformations

- Lines 37: computes `hz_per_mt` using `hz_per_mt=abs(spin('E'))/(2*pi)*1e-3`.
- Lines 40: computes `centre` using `centre=lower(parameters.centre)`.
- Lines 41: computes `electron` using `electron='E'`.
- Lines 46: computes `gvals` using `gvals=[2.3463 1.8438 1.7045]`.
- Lines 47: computes `Aco` using `Aco=[8.86 6.43 5.82]`.
- Lines 48: computes `alpha` using `alpha=29`.
- Lines 58: computes `yaxis` using `yaxis=[0;1;1]/sqrt(2)`.
- Lines 59: computes `xbase` using `xbase=[1;0;0]`.
- Lines 60: computes `zbase` using `zbase=cross(xbase,yaxis)`.
- Lines 61: computes `xaxis` using `xaxis=cosd(alpha)*xbase+sind(alpha)*zbase`.
- Lines 62: computes `zaxis` using `zaxis=cross(xaxis,yaxis)`.
- Lines 63: computes `frame` using `frame=[xaxis(:) yaxis(:) zaxis(:)]`.
- Lines 66: computes `[frame,~]` using `[frame,~]=qr(frame,0)`.
- Lines 70: computes `frame(:,3)` using `frame(:,3)=-frame(:,3)`.
- Lines 74: computes `gmat` using `gmat=((frame)*diag(gvals)*(frame)')`.
- Lines 75: computes `nucleus` using `nucleus=struct('iso','59Co','A',((frame)*diag(Aco*hz_per_mt)*(frame)'))`.
- Lines 80: computes `C` using `C=rotmat_align([1 1 1],[0 0 1])`.
- Lines 88: computes `sys.isotopes` using `sys.isotopes={electron,nucleus.iso}`.

### Local helper functions

- Line 97: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if(~isstruct(parameters))`.
  - Representative operation: `error('parameters must be a structure.')`.

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'o4' or 'nlo2'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Cobalt-related defect spin system for diamond. Syntax:
- [sys,inter]=diamond_co(parameters)
- Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),
- parameters is a structure with the following fields:
- .centre -'o4' or 'nlo2'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure
- Check input count
- Check consistency
- Set field-unit conversion constants

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `spin()`, `lower()`, `cross()`, `cosd()`, `sind()`, `xaxis()`, `yaxis()`, `zaxis()`, `frame()`, `rotmat_align()`, `isstruct()`, `isfield()`, `ischar()`, `any()`, `strcmpi()`.
