# kernel/conventions/transforms/eeqq2nqi.m

- Signature: `Q=eeqq2nqi(C_q,eta_q,I,eulers)`

## Purpose

Converts the C_q and eta_q quadrupolar interaction specification convention into a 3x3 interaction matrix in Hz. Syntax: Q=eeqq2nqi(C_q,eta_q,I,euler_angles)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.
- Quadrupolar physics is relevant: nuclei with spin > 1/2 interact with the electric field gradient tensor, introducing second-rank anisotropy, asymmetry, and overtone or MQ phenomena.

## Numerical / algorithmic content

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
