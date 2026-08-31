# kernel/conventions/transforms/euler2dcm.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/euler2dcm.m`
- Signature: `R=euler2dcm(arg1,arg2,arg3)`
- Total lines: 79

## Purpose

Converts Euler angles (ZYZ active convention) into a direction cosine matrix. Syntax: R=euler2dcm(alpha,beta,gamma) OR R=euler2dcm([alpha beta gamma])

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Adapt to the input style; implemented by `if nargin==1`.
- Lines 32-33: Assume that a single input is a 3-vector; implemented by `alp=arg1(1); bet=arg1(2); gam=arg1(3)`.
- Lines 35-36: Assume that three inputs are Euler angles; implemented by `alp=arg1; bet=arg2; gam=arg3`.
- Lines 38-39: Bomb out in all other cases; implemented by `error('incorrect number of input arguments.')`.
- Lines 42-43: Check consistency; implemented by `grumble(alp,bet,gam)`.
- Lines 45-47: Build the individual rotation matrices, as per Brink & Satchler, Fig 1a; implemented by `R_alpha=[cos(alp) -sin(alp) 0`.
- Lines 57-58: Build the direction cosine matrix; implemented by `R=R_alpha*R_beta*R_gamma`.

### Control flow inferred from the code

- Line 31: conditional branch on `nargin==1`.

### Key state/data transformations

- Lines 33: computes `alp` using `alp=arg1(1); bet=arg1(2); gam=arg1(3)`.
- Lines 47: computes `R_alpha` using `R_alpha=[cos(alp) -sin(alp) 0`.
- Lines 50: computes `R_beta` using `R_beta= [cos(bet) 0 sin(bet)`.
- Lines 53: computes `R_gamma` using `R_gamma=[cos(gam) -sin(gam) 0`.
- Lines 58: computes `R` using `R=R_alpha*R_beta*R_gamma`.

### Local helper functions

- Line 63: `grumble()` — `function grumble(alp,bet,gam)`.
  - Representative operation: `if (~isnumeric(alp))||(~isnumeric(bet))||(~isnumeric(gam))`.
  - Representative operation: `error('all inputs must be numeric.')`.

## Parameters / inputs

- alpha,beta,gamma -Euler angles in radians (ZYZ
- active convention)

## Outputs

- R -direction cosine matrix
- Note: the resulting rotation matrix is to be used as follows:
- v=R*v (for 3x1 vectors)
- A=R*A*R' (for 3x3 interaction tensors)

## Implementation structure

- Converts Euler angles (ZYZ active convention) into a direction
- cosine matrix. Syntax:
- R=euler2dcm(alpha,beta,gamma)
- R=euler2dcm([alpha beta gamma])
- alpha,beta,gamma - Euler angles in radians (ZYZ
- active convention)
- R - direction cosine matrix
- Note: the resulting rotation matrix is to be used as follows:
- v=R*v (for 3x1 vectors)
- A=R*A*R' (for 3x3 interaction tensors)
- Adapt to the input style
- Assume that a single input is a 3-vector

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `arg1()`, `grumble()`.
