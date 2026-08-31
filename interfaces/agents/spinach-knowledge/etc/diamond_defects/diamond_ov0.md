# etc/diamond_defects/diamond_ov0.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/diamond_defects/diamond_ov0.m`
- Signature: `[sys,inter]=diamond_ov0(parameters)`
- Total lines: 155

## Purpose

Neutral oxygen-vacancy (OV0, WAR5) centre ground state spin system for diamond. Syntax: [sys,inter]=diamond_ov0(parameters) Magnetic parameters from Tables 9-2 and 9-3 of: B.L. Cann, Magnetic Resonance Studies of Point Defects in Diamond, PhD thesis, University of Warwick (2009), The zero-field splitting is confirmed at 4 K, and the defect is assigned to OV0, in: S. Mukherjee et al., Phys. Rev. B 114, 074105 (2026), 

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 52-53: Check input count; implemented by `if nargin~=1`.
- Lines 57-58: Check consistency; implemented by `grumble(parameters)`.
- Lines 60-63: Set the trigonal principal-axis frame; implemented by `frame=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 65-66: Rotation matrix for orientation; implemented by `switch parameters.orientation`.
- Lines 82-83: Complain and bomb out; implemented by `error('unknown orientation specification.')`.
- Lines 87-88: Define spin system isotopes and labels; implemented by `sys.isotopes={'E3','13C','13C','13C','13C'}`.
- Lines 93-94: Electron g-tensor; implemented by `inter.zeeman.matrix{1}=C*frame*diag([2.0025 2.0025 2.0029])*frame'*C'`.
- Lines 96-97: Electron ZFS tensor; implemented by `inter.coupling.matrix{1,1}=C*frame*zfs2mat(2888e6,0,0,0,0)*frame'*C'`.
- Lines 99-100: 13C hyperfine eigenvalues, MHz; implemented by `hfc_vals=[197.4 117.3 118.2`.
- Lines 105-106: 13C hyperfine principal-axis polar angles from [001], degrees; implemented by `hfc_theta=[53.86 143.86 90.00`.
- Lines 111-112: 13C hyperfine principal-axis azimuths from [100], degrees; implemented by `hfc_phi=[225.26 225.00 135.00`.
- Lines 117-118: Build and rotate the 13C hyperfine tensors; implemented by `for n=1:size(hfc_vals,1)`.
- Lines 120-123: Convert principal-axis directions from polar angles; implemented by `axis_xyz=[sind(hfc_theta(n,:)).*cosd(hfc_phi(n,:)); sind(hfc_theta(n,:)).*sind(hfc_phi(n,:)); cosd(hfc_theta(n,:))]`.
- Lines 125-126: Build the crystal-frame hyperfine tensor; implemented by `hfc_cub=axis_xyz*diag(1e6*hfc_vals(n,:))*axis_xyz'`.
- Lines 128-129: Symmetrise and rotate into the requested orientation; implemented by `inter.coupling.matrix{1,n+1}=C*((hfc_cub+hfc_cub')/2)*C'`.

### Control flow inferred from the code

- Line 53: conditional branch on `nargin~=1`.
- Line 66: dispatches on `parameters.orientation`; cases `'100'`, `'110'`, `'111'`.
- Line 118: `for` loop over `n=1:size(hfc_vals,1)`.

### Key state/data transformations

- Lines 61-63: computes `frame` using `frame=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 1/sqrt(2) -1/sqrt(6) 1/sqrt(3); 0 2/sqrt(6) 1/sqrt(3)]`.
- Lines 70: computes `C` using `C=rotmat_align([1 0 0],[0 0 1])`.
- Lines 88: computes `sys.isotopes` using `sys.isotopes={'E3','13C','13C','13C','13C'}`.
- Lines 89: computes `sys.labels` using `sys.labels={'OV0','C_a','C_g','C_l','C_c'}`.
- Lines 90: computes `inter.zeeman.matrix` using `inter.zeeman.matrix=cell(1,numel(sys.isotopes))`.
- Lines 91: computes `inter.coupling.matrix` using `inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes))`.
- Lines 94: computes `inter.zeeman.matrix{1}` using `inter.zeeman.matrix{1}=C*frame*diag([2.0025 2.0025 2.0029])*frame'*C'`.
- Lines 97: computes `inter.coupling.matrix{1,1}` using `inter.coupling.matrix{1,1}=C*frame*zfs2mat(2888e6,0,0,0,0)*frame'*C'`.
- Lines 100: computes `hfc_vals` using `hfc_vals=[197.4 117.3 118.2`.
- Lines 106: computes `hfc_theta` using `hfc_theta=[53.86 143.86 90.00`.
- Lines 112: computes `hfc_phi` using `hfc_phi=[225.26 225.00 135.00`.
- Lines 121-123: computes `axis_xyz` using `axis_xyz=[sind(hfc_theta(n,:)).*cosd(hfc_phi(n,:)); sind(hfc_theta(n,:)).*sind(hfc_phi(n,:)); cosd(hfc_theta(n,:))]`.
- Lines 126: computes `hfc_cub` using `hfc_cub=axis_xyz*diag(1e6*hfc_vals(n,:))*axis_xyz'`.
- Lines 129: computes `inter.coupling.matrix{1,n+1}` using `inter.coupling.matrix{1,n+1}=C*((hfc_cub+hfc_cub')/2)*C'`.

### Local helper functions

- Line 136: `grumble()` — `function grumble(parameters)`.
  - Representative operation: `if(~isstruct(parameters))`.
  - Representative operation: `error('parameters must be a structure.')`.

## Parameters / inputs

- a structure (parameters.*) with the following field:
- .orientation -'111', '110', or '100' crystal plane normal
- aligned with the magnetic field

## Outputs

- sys -Spinach system specification structure
- inter -Spinach interaction specification structure

## Implementation structure

- Neutral oxygen-vacancy (OV0, WAR5) centre ground state spin system
- for diamond. Syntax:
- [sys,inter]=diamond_ov0(parameters)
- Magnetic parameters from Tables 9-2 and 9-3 of:
- B.L. Cann, Magnetic Resonance Studies of Point Defects in
- Diamond, PhD thesis, University of Warwick (2009),
- The zero-field splitting is confirmed at 4 K, and the defect is
- assigned to OV0, in:
- S. Mukherjee et al., Phys. Rev. B 114, 074105 (2026),
- The centre has S=1 and C3v symmetry; the electron Zeeman and zero-
- field splitting tensors are axial about the trigonal axis, and the
- rhombicity is zero within the experimental error. Oxygen is left out

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `rotmat_align()`, `zfs2mat()`, `sind()`, `hfc_theta()`, `cosd()`, `hfc_phi()`, `hfc_vals()`, `isstruct()`, `isfield()`, `ischar()`, `ismember()`.
