# kernel/conventions/transforms/eeqq2nqi.m

- Signature: `Q=eeqq2nqi(C_q,eta_q,I,eulers)`

## Purpose

Converts the C_q and eta_q quadrupolar interaction specification convention into a 3x3 interaction matrix in Hz. Syntax: Q=eeqq2nqi(C_q,eta_q,I,euler_angles)

## Physical / mathematical content
The function constructs a quadrupolar coupling tensor from the coupling constant C_q, asymmetry parameter eta_q, and spin quantum number I. Its principal-axis values are XX = -C_q(1-eta_q)/[4I(2I-1)], YY = -C_q(1+eta_q)/[4I(2I-1)], and ZZ = C_q/[2I(2I-1)]. The Euler angles orient that principal-axis tensor relative to the lab frame.

## Numerical / algorithmic content
The function obtains a rotation matrix R from euler2dcm(eulers), then computes Q = R*diag([XX YY ZZ])*R'. It removes any residual isotropic component by subtracting eye(3)*trace(Q)/3 and symmetrizes the result as (Q+Q')/2.

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
The function first calls a local grumble routine to validate the inputs: all must be numeric and real; eulers must have three elements; C_q, eta_q, and I must each be scalar; and I must be an integer or half-integer at least 1. It then calculates the principal-axis values, rotates the diagonal tensor, and cleans up numerical trace and symmetry errors before returning Q.
