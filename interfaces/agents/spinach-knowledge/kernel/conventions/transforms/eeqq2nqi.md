# kernel/conventions/transforms/eeqq2nqi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/eeqq2nqi.m`
- Signature: `Q=eeqq2nqi(C_q,eta_q,I,eulers)`
- Total lines: 76

## Purpose

Converts the C_q and eta_q quadrupolar interaction specification convention into a 3x3 interaction matrix in Hz. Syntax: Q=eeqq2nqi(C_q,eta_q,I,euler_angles)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(C_q,eta_q,I,eulers)`.
- Lines 38-39: Get the eigenvalues; implemented by `XX=-C_q*(1-eta_q)/(4*I*(2*I-1))`.
- Lines 43-44: Get the rotation matrix; implemented by `R=euler2dcm(eulers)`.
- Lines 46-47: Get the quadrupole tensor matrix; implemented by `Q=R*diag([XX YY ZZ])*R'`.
- Lines 49-50: Clean up small rounding errors; implemented by `Q=Q-eye(3)*trace(Q)/3; Q=(Q+Q')/2`.

### Key state/data transformations

- Lines 39: computes `XX` using `XX=-C_q*(1-eta_q)/(4*I*(2*I-1))`.
- Lines 40: computes `YY` using `YY=-C_q*(1+eta_q)/(4*I*(2*I-1))`.
- Lines 41: computes `ZZ` using `ZZ=+C_q/(2*I*(2*I-1))`.
- Lines 44: computes `R` using `R=euler2dcm(eulers)`.
- Lines 47: computes `Q` using `Q=R*diag([XX YY ZZ])*R'`.

### Local helper functions

- Line 55: `grumble()` — `function grumble(C_q,eta_q,I,eulers)`.
  - Representative operation: `if (~isnumeric(C_q))||(~isnumeric(eta_q))||(~isnumeric(I))||(~isnumeric(eulers))`.
  - Representative operation: `error('all inputs must be numeric.')`.

## Parameters / inputs

- C_q -quadrupolar coupling constant e^2*q*Q/h
- in Hz
- eta_q -quadrupolar tensor asymmetry parameter
- I -spin quantum number
- euler_angles -vector of three Euler angles in radians,
- giving the orientation of the principal
- axis frame relative to the lab frame.

## Outputs

- Q -quadrupolar coupling tensor as a 3x3
- matrix in Hz
- Note: the denominator contains spin quantum number squared,
- meaning that the actual 3x3 interaction tensor falls
- off sharply with the nuclear spin quantum number for
- the same anisotropy parameter.

## Implementation structure

- Converts the C_q and eta_q quadrupolar interaction specification
- convention into a 3x3 interaction matrix in Hz. Syntax:
- Q=eeqq2nqi(C_q,eta_q,I,euler_angles)
- C_q -quadrupolar coupling constant e^2*q*Q/h
- in Hz
- eta_q -quadrupolar tensor asymmetry parameter
- I -spin quantum number
- euler_angles -vector of three Euler angles in radians,
- giving the orientation of the principal
- axis frame relative to the lab frame.
- Q -quadrupolar coupling tensor as a 3x3
- matrix in Hz

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `euler2dcm()`.
