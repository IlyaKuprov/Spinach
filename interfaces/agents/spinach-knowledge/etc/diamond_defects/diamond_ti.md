# etc/diamond_defects/diamond_ti.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_ti.m`
- Signature: `[sys,inter]=diamond_ti(parameters)`
- Total lines: 186

## Purpose

Titanium-related defect spin system for diamond. Syntax: [sys,inter]=diamond_ti(parameters) Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `diamond_frame_xyz()`, `diamond_frame_alpha()`, `diamond_frame_xz()`, `round()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check input count; implemented by `if nargin~=1`.
- Lines 36-37: Check consistency; implemented by `grumble(parameters)`.
- Lines 39-40: Set irrelevant carbon count to zero; implemented by `if ~strcmpi(parameters.centre,'ok1')`.
- Lines 44-45: Set field-unit conversion constants; implemented by `hz_per_mt=abs(spin('E'))/(2*pi)*1e-3`.
- Lines 47-48: Select the titanium centre; implemented by `centre=lower(parameters.centre)`.
- Lines 51-52: Select centre-specific magnetic parameters; implemented by `switch centre`.
- Lines 70-71: Add the nitrogen hyperfine tensor; implemented by `nuclei{end+1}=struct('iso','14N','A',((aframe)*diag(An*hz_per_mt)*(aframe)'))`.
- Lines 73-74: Add the titanium isotope if requested; implemented by `titanium=parameters.titanium`.
- Lines 79-80: Add reported OK1 nearest-neighbour carbons; implemented by `if strcmp(centre,'ok1')&&parameters.n_13c>0`.
- Lines 89-90: Build the Spinach structures; implemented by `switch parameters.orientation`.

### Control flow inferred from the code

- Line 32: conditional branch on `nargin~=1`.
- Line 40: conditional branch on `~strcmpi(parameters.centre,'ok1')`.
- Line 52: dispatches on `centre`; cases `'n3'`, `'ok1'`.
- Line 75: conditional branch on `~strcmp(titanium,'none')`.
- Line 80: conditional branch on `strcmp(centre,'ok1')&&parameters.n_13c>0`.
- Line 83: `for` loop over `n=1:parameters.n_13c`.
- Line 90: dispatches on `parameters.orientation`; cases `'111'`, `'110'`, `'100'`.
- Line 101: `for` loop over `n=1:numel(nuclei)`.
- Line 107: `for` loop over `n=1:numel(nuclei)`.
- Line 108: conditional branch on `isfield(nuclei{n},'A')&&norm(nuclei{n}.A,2)>0`.

### Key state/data transformations

- Lines 41: computes `parameters.n_13c` using `parameters.n_13c=0`.
- Lines 45: computes `hz_per_mt` using `hz_per_mt=abs(spin('E'))/(2*pi)*1e-3`.
- Lines 48: computes `centre` using `centre=lower(parameters.centre)`.
- Lines 49: computes `electron` using `electron='E'; nuclei={}`.
- Lines 54: computes `gvals` using `gvals=[2.0022 2.0025 2.0020]`.
- Lines 55: computes `An` using `An=[0.11 0.15 0.11]`.
- Lines 56: computes `Ati` using `Ati=[0.28 0.40 0.28]`.
- Lines 57: computes `g_alpha` using `g_alpha=32; A_alpha=26`.
- Lines 66: computes `gframe` using `gframe=diamond_frame_alpha(g_alpha)`.
- Lines 67: computes `aframe` using `aframe=diamond_frame_alpha(A_alpha)`.
- Lines 68: computes `gmat` using `gmat=((gframe)*diag(gvals)*(gframe)')`.
- Lines 71: computes `nuclei{end+1}` using `nuclei{end+1}=struct('iso','14N','A',((aframe)*diag(An*hz_per_mt)*(aframe)'))`.
- Lines 74: computes `titanium` using `titanium=parameters.titanium`.
- Lines 81: computes `cframe{1}` using `cframe{1}=diamond_frame_xz([1 1 0],[1 -1 -1])`.
- Lines 82: computes `cframe{2}` using `cframe{2}=diamond_frame_xz([1 1 0],[-1 1 -1])`.
- Lines 84: computes `Cmat` using `Cmat=((cframe{n})*diag([2.62 2.62 4.38]*hz_per_mt)*(cframe{n})')`.
- Lines 92: computes `C` using `C=rotmat_align([1 1 1],[0 0 1])`.
- Lines 100: computes `sys.isotopes` using `sys.isotopes={electron}`.

### Local helper functions

- Line 116: `diamond_frame_xyz()` — `function frame=diamond_frame_xyz(xaxis,yaxis,zaxis)`. Frame used in the Nadolinny nickel and titanium tables
  - Representative operation: `frame=[xaxis(:) yaxis(:) zaxis(:)]`.
  - Representative operation: `[frame,~]=qr(frame,0)`.
- Line 125: `diamond_frame_alpha()` — `function frame=diamond_frame_alpha(alpha)`. Principal-axis frame from a specified x axis and z axis
  - Representative operation: `xaxis=[1;-1;0]/sqrt(2)`.
  - Representative operation: `ybase=[1;1;0]/sqrt(2)`.
- Line 135: `diamond_frame_xz()` — `function frame=diamond_frame_xz(xaxis,zaxis)`. Consistency enforcement
  - Representative operation: `xaxis=xaxis(:)/norm(xaxis,2)`.
  - Representative operation: `zaxis=zaxis(:)-xaxis*dot(xaxis,zaxis(:))`.
- Line 144: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if ~isstruct(parameters)`.
  - Representative operation: `error('parameters must be a structure.')`.

## Parameters / inputs

- parameters is a structure with the following fields:
- .centre -'n3' or 'ok1'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .titanium -titanium isotope label, or 'none'
- .n_13c -number of reported 13C hyperfine couplings
- to include, from 0 to 2; applies only to OK1

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Titanium-related defect spin system for diamond. Syntax:
- [sys,inter]=diamond_ti(parameters)
- Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),
- parameters is a structure with the following fields:
- .centre -'n3' or 'ok1'
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field
- .titanium -titanium isotope label, or 'none'
- .n_13c -number of reported 13C hyperfine couplings
- to include, from 0 to 2; applies only to OK1
- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `strcmpi()`, `spin()`, `lower()`, `diamond_frame_alpha()`, `strcmp()`, `diamond_frame_xz()`, `rotmat_align()`, `isfield()`, `diamond_frame_xyz()`, `xaxis()`, `yaxis()`, `zaxis()`, `frame()`, `cosd()`, `sind()`.
