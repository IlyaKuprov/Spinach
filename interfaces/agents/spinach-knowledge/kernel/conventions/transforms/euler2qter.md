# kernel/conventions/transforms/euler2qter.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/euler2qter.m`
- Signature: `q=euler2qter(arg1,arg2,arg3)`
- Total lines: 73

## Purpose

Converts Euler angles (ZYZ active convention) into a unit quaternion in the active convention, matching euler2dcm.m function. Syntax: q=euler2qter(alpha,beta,gamma) OR q=euler2qter([alpha beta gamma])

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 33-34: Adapt to the input style; implemented by `if nargin==1`.
- Lines 36-37: Assume that a single input is a 3-vector; implemented by `alp=arg1(1); bet=arg1(2); gam=arg1(3)`.
- Lines 41-42: Assume that three inputs are Euler angles; implemented by `alp=arg1; bet=arg2; gam=arg3`.
- Lines 46-47: Bomb out in all other cases; implemented by `error('incorrect number of input arguments.')`.
- Lines 51-52: Check consistency; implemented by `grumble(alp,bet,gam)`.
- Lines 54-55: Compute the quaternion; implemented by `q.u=cos(bet/2).*cos((alp+gam)/2)`.

### Control flow inferred from the code

- Line 34: conditional branch on `nargin==1`.

### Key state/data transformations

- Lines 37: computes `alp` using `alp=arg1(1); bet=arg1(2); gam=arg1(3)`.
- Lines 55: computes `q.u` using `q.u=cos(bet/2).*cos((alp+gam)/2)`.
- Lines 56: computes `q.i` using `q.i=sin(bet/2).*sin((gam-alp)/2)`.
- Lines 57: computes `q.j` using `q.j=sin(bet/2).*cos((gam-alp)/2)`.
- Lines 58: computes `q.k` using `q.k=cos(bet/2).*sin((alp+gam)/2)`.

### Local helper functions

- Line 63: `grumble()` — `function grumble(alp,bet,gam)`.
  - Representative operation: `if (~isnumeric(alp))||(~isreal(alp))||(~iscolumn(alp))|| (~isnumeric(bet))||(~isreal(bet))||(~iscolumn(bet))|| (~isnumeric(gam))||(~isreal(gam))||(~iscolumn(gam))`.
  - Representative operation: `(~isnumeric(bet))||(~isreal(bet))||(~iscolumn(bet))|| (~isnumeric(gam))||(~isreal(gam))||(~iscolumn(gam))`.

## Parameters / inputs

- alpha,beta,gamma -Euler angles in radians (ZYZ active
- convention), scalars or column vec-
- tors of equal length

## Outputs

- q -structure with four fields q.u, q.i, q.j, q.k giving
- the four components of the quaternion; for column
- vector inputs each field is a column vector
- Note: the quaternion returned represents the same rotation
- as euler2dcm(alpha,beta,gamma); it is converted into
- that matrix by qter2dcm.m function.

## Implementation structure

- Converts Euler angles (ZYZ active convention) into a unit
- quaternion in the active convention, matching euler2dcm.m
- function. Syntax:
- q=euler2qter(alpha,beta,gamma)
- q=euler2qter([alpha beta gamma])
- alpha,beta,gamma -Euler angles in radians (ZYZ active
- convention), scalars or column vec-
- tors of equal length
- q -structure with four fields q.u, q.i, q.j, q.k giving
- the four components of the quaternion; for column
- vector inputs each field is a column vector
- Note: the quaternion returned represents the same rotation

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `arg1()`, `grumble()`, `iscolumn()`.
