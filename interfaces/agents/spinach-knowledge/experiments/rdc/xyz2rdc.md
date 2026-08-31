# experiments/rdc/xyz2rdc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/experiments/rdc/xyz2rdc.m`
- Signature: `rdc=xyz2rdc(spin_a,spin_b,xyz_a,xyz_b,order_spec)`
- Total lines: 86

## Purpose

Converts Cartesian coordinates of a pair of nuclei and an order matrix into residual dipolar coupling; only hetero- nuclear spin pairs are supported. Syntax: rdc=xyz2rdc(spin_a,spin_b,xyz_a,xyz_b,chi)

## Physical / mathematical content

- Residual-dipolar-coupling experiment and analysis routines. These files use partial ordering, Saupe tensors, and molecular-frame geometry to connect internuclear vectors with observed couplings.
- The spin physics includes through-space magnetic dipole-dipole coupling, a rank-2 anisotropic interaction with strong orientation dependence and characteristic secular/non-secular structure.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 31-32: Check consistency; implemented by `grumble(spin_a,spin_b,xyz_a,xyz_b,order_spec)`.
- Lines 34-35: Choose the path; implemented by `switch order_spec{2}`.
- Lines 37-38: Saupe order matrix; implemented by `case 'saupe'`.
- Lines 40-41: Get the dipole coupling tensor in rad/s; implemented by `[~,~,~,~,D]=xyz2dd(xyz_a,xyz_b,spin_a,spin_b)`.
- Lines 43-44: Get weak heteronuclear RDC in Hz; implemented by `rdc=(2/3)*trace(order_spec{1}*D)/(2*pi)`.
- Lines 46-47: For future options; implemented by `otherwise`.
- Lines 49-50: Complain and bomb out; implemented by `error('unknown order specification type')`.

### Control flow inferred from the code

- Line 35: dispatches on `order_spec{2}`; cases `'saupe'`.

### Key state/data transformations

- Lines 41: computes `[~,~,~,~,D]` using `[~,~,~,~,D]=xyz2dd(xyz_a,xyz_b,spin_a,spin_b)`.
- Lines 44: computes `rdc` using `rdc=(2/3)*trace(order_spec{1}*D)/(2*pi)`.

### Local helper functions

- Line 57: `grumble()` — `function grumble(spin_a,spin_b,xyz_a,xyz_b,order_spec)`.
  - Representative operation: `if (~iscell(order_spec))||(~ismember(numel(order_spec),[2 4]))`.
  - Representative operation: `error('order_spec must be a cell array with 2 or 4 elements.')`.

## Parameters / inputs

- spin_a, spin_b -character strings indicating spin
- type, for example '13C'
- xyz_a, xyz_b -three-element vectors specifying
- Cartesian coordinates of the two
- spins in Angstroms
- order_spec -{S,'saupe'} uses Saupe order mat-
- rix, S is a traceless symmetric
- 3x3 matrix, dimensionless

## Outputs

- rdc -residual dipolar coupling in the
- heteronuclear case, Hz

## Implementation structure

- Converts Cartesian coordinates of a pair of nuclei and an
- order matrix into residual dipolar coupling; only hetero-
- nuclear spin pairs are supported. Syntax:
- rdc=xyz2rdc(spin_a,spin_b,xyz_a,xyz_b,chi)
- spin_a, spin_b -character strings indicating spin
- type, for example '13C'
- xyz_a, xyz_b -three-element vectors specifying
- Cartesian coordinates of the two
- spins in Angstroms
- order_spec -{S,'saupe'} uses Saupe order mat-
- rix, S is a traceless symmetric
- 3x3 matrix, dimensionless

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `xyz2dd()`, `iscell()`, `ismember()`, `ismatrix()`, `any()`, `eps()`, `ischar()`, `strcmp()`.
