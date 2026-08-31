# kernel/conventions/transforms/euler_sup.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/euler_sup.m`
- Signature: `rot_cmp=euler_sup(rot_one,rot_two)`
- Total lines: 107

## Purpose

Superposition of ZYZ active Euler rotations. Syntax: rot_cmp=euler_sup(rot_one,rot_two)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(rot_one,rot_two)`.
- Lines 35-36: Cast to row vectors; implemented by `rot_one=rot_one(:)'`.
- Lines 39-40: Wrap both angle sets; implemented by `rot_one=wrapToPi(rot_one)`.
- Lines 43-44: Build DCMs for both; implemented by `dcm_one=euler2dcm(rot_one)`.
- Lines 47-48: Compose the rotations; implemented by `dcm_comp=dcm_two*dcm_one`.
- Lines 50-51: Handle identity rotation special case; implemented by `if all(abs(diag(dcm_comp)-[1;1;1])<1e-12)`.
- Lines 53-54: Identity rotation; implemented by `rot_cmp=[pi/2 0 -pi/2]`.
- Lines 56-59: Handle the singular beta=pi branch; implemented by `elseif (abs(dcm_comp(3,3)+1)<1e-12)&& (abs(rot_one(1)-rot_two(1))<1e-12)&& (abs(rot_one(3)-rot_two(3))<1e-12)`.
- Lines 61-62: Clamp the cosine argument; implemented by `dcm_11=max(min(dcm_comp(1,1),1),-1)`.
- Lines 64-65: Recover the third angle; implemented by `gam=(acos(dcm_11)-pi)/2`.
- Lines 67-68: Return special case; implemented by `rot_cmp=[-gam pi gam]`.
- Lines 72-73: Process the rest normally; implemented by `rot_cmp=dcm2euler(dcm_comp)`.
- Lines 77-78: Wrap output into (-pi,pi]; implemented by `rot_cmp=wrapToPi(rot_cmp)`.

### Control flow inferred from the code

- Line 51: conditional branch on `all(abs(diag(dcm_comp)-[1;1;1])<1e-12)`.

### Key state/data transformations

- Lines 36: computes `rot_one` using `rot_one=rot_one(:)'`.
- Lines 37: computes `rot_two` using `rot_two=rot_two(:)'`.
- Lines 44: computes `dcm_one` using `dcm_one=euler2dcm(rot_one)`.
- Lines 45: computes `dcm_two` using `dcm_two=euler2dcm(rot_two)`.
- Lines 48: computes `dcm_comp` using `dcm_comp=dcm_two*dcm_one`.
- Lines 54: computes `rot_cmp` using `rot_cmp=[pi/2 0 -pi/2]`.
- Lines 62: computes `dcm_11` using `dcm_11=max(min(dcm_comp(1,1),1),-1)`.
- Lines 65: computes `gam` using `gam=(acos(dcm_11)-pi)/2`.

### Local helper functions

- Line 83: `grumble()` — `function grumble(rot_one,rot_two)`. "My mathematics work is proceeding beyond my wildest
  - Representative operation: `if (~isnumeric(rot_one))||(~isreal(rot_one))||(~isvector(rot_one))|| (numel(rot_one)~=3)||any(~isfinite(rot_one(:)))`.
  - Representative operation: `(numel(rot_one)~=3)||any(~isfinite(rot_one(:)))`.

## Parameters / inputs

- rot_one -first Euler angle set [alpha beta gamma],
- radians, ZYZ active convention
- rot_two -second Euler angle set [alpha beta gamma],
- radians, ZYZ active convention

## Outputs

- rot_cmp -row vector [alpha beta gamma] of the
- composite rotation, radians
- Note: rotations are applied in the supplied order
- v_rot=R_two*R_one*v
- therefore the composite matrix is
- R_comp=R_two*R_one

## Implementation structure

- Superposition of ZYZ active Euler rotations. Syntax:
- rot_cmp=euler_sup(rot_one,rot_two)
- rot_one -first Euler angle set [alpha beta gamma],
- radians, ZYZ active convention
- rot_two -second Euler angle set [alpha beta gamma],
- rot_cmp -row vector [alpha beta gamma] of the
- composite rotation, radians
- Note: rotations are applied in the supplied order
- v_rot=R_two*R_one*v
- therefore the composite matrix is
- R_comp=R_two*R_one
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `rot_one()`, `rot_two()`, `wrapToPi()`, `euler2dcm()`, `all()`, `elseif()`, `dcm_comp()`, `acos()`, `dcm2euler()`, `isvector()`, `any()`.
