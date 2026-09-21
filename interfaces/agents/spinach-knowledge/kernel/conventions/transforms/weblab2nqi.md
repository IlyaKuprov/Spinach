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
