# kernel/conventions/transforms/weblab2nqi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/weblab2nqi.m`
- Signature: `varargout=weblab2nqi(C_q,eta_q,I,alpha,theta,phi)`
- Total lines: 118

## Purpose

Converts the Weblab one-cone model parameters (see weblab_cone.png) into NQI tensors used by Spinach. Syntax: [Q1,Q2]=weblab2nqi(C_q,eta_q,I,alpha,theta,phi) [Q1,Q2,Q3]=weblab2nqi(C_q,eta_q,I,alpha,theta,phi) [Q1,Q2,Q3,Q4]=weblab2nqi(C_q,eta_q,I,alpha,theta) [Q1,Q2,Q3,Q4,Q5,Q6]=weblab2nqi(C_q,eta_q,I,alpha,theta)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 39-40: Four-and six-output modes take no phi; implemented by `if nargin<6, phi=[]; end`.
- Lines 42-43: Check consistency; implemented by `grumble(C_q,eta_q,I,alpha,theta,phi,nargin,nargout)`.
- Lines 45-46: Translate conventions and call eeqq2nqi; implemented by `switch nargout`.
- Lines 50-51: Two output arguments; implemented by `varargout{1}=eeqq2nqi(C_q,eta_q,I,[-phi/2 theta alpha])`.
- Lines 56-57: Three output arguments; implemented by `varargout{1}=eeqq2nqi(C_q,eta_q,I,[-phi theta alpha])`.
- Lines 63-64: Four output arguments, fixed azimuth grid; implemented by `varargout{1}=eeqq2nqi(C_q,eta_q,I,[0*pi/2 theta alpha])`.
- Lines 71-72: Six output arguments, fixed azimuth grid; implemented by `varargout{1}=eeqq2nqi(C_q,eta_q,I,[0*pi/3 theta alpha])`.
- Lines 81-82: Complain and bomb out; implemented by `error('incorrect number of output arguments')`.

### Control flow inferred from the code

- Line 40: conditional branch on `nargin<6, phi=[]; end`.
- Line 46: dispatches on `nargout`; cases `2`, `3`, `4`, `6`.

### Key state/data transformations

- Lines 51: computes `varargout{1}` using `varargout{1}=eeqq2nqi(C_q,eta_q,I,[-phi/2 theta alpha])`.
- Lines 52: computes `varargout{2}` using `varargout{2}=eeqq2nqi(C_q,eta_q,I,[+phi/2 theta alpha])`.
- Lines 59: computes `varargout{3}` using `varargout{3}=eeqq2nqi(C_q,eta_q,I,[+phi theta alpha])`.
- Lines 67: computes `varargout{4}` using `varargout{4}=eeqq2nqi(C_q,eta_q,I,[3*pi/2 theta alpha])`.
- Lines 76: computes `varargout{5}` using `varargout{5}=eeqq2nqi(C_q,eta_q,I,[4*pi/3 theta alpha])`.
- Lines 77: computes `varargout{6}` using `varargout{6}=eeqq2nqi(C_q,eta_q,I,[5*pi/3 theta alpha])`.

### Local helper functions

- Line 89: `grumble()` — `function grumble(C_q,eta_q,I,alpha,theta,phi,n_ins,n_outs)`.
  - Representative operation: `if ((n_outs==2)||(n_outs==3))&&(n_ins~=6)`.
  - Representative operation: `error('phi is required in the two- and three-output modes.')`.

## Parameters / inputs

- C_q -quadrupolar coupling constant e^2*q*Q/h
- in Hz
- eta_q -quadrupolar tensor asymmetry parameter
- I -spin quantum number
- alpha
- theta
- phi -the three angles of Weblab cone model
- (see weblab_cone.png), in radians; the
- four-and six-site modes place the sites
- on fixed azimuth grids, [0 1 2 3]*pi/2
- and [0 1 2 3 4 5]*pi/3 respectively, and
- must therefore be called without phi

## Outputs

- Q1,Q2,... -quadrupolar coupling tensors for the two,
- three, four, or six sites as 3x3 matrices
- in Hz

## Implementation structure

- Converts the Weblab one-cone model parameters (see weblab_cone.png)
- into NQI tensors used by Spinach. Syntax:
- [Q1,Q2]=weblab2nqi(C_q,eta_q,I,alpha,theta,phi)
- [Q1,Q2,Q3]=weblab2nqi(C_q,eta_q,I,alpha,theta,phi)
- [Q1,Q2,Q3,Q4]=weblab2nqi(C_q,eta_q,I,alpha,theta)
- [Q1,Q2,Q3,Q4,Q5,Q6]=weblab2nqi(C_q,eta_q,I,alpha,theta)
- C_q -quadrupolar coupling constant e^2*q*Q/h
- in Hz
- eta_q -quadrupolar tensor asymmetry parameter
- I -spin quantum number
- alpha
- theta

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `eeqq2nqi()`, `isscalar()`.
