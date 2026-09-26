# kernel/conventions/transforms/castep2nqi.m

- Signature: `nqi=castep2nqi(V,Q,I)`

## Purpose

Converts CASTEP EFG tensor (it is printed in atomic units) to NQI 3x3 tensor in Hz that is required by Spinach. Syntax: nqi=castep2nqi(V,Q,I)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- V -EFG tensor from CASTEP output, a.u.
- Q -nuclear quadrupole moment, barn
- I -nuclear spin quantum number

## Outputs

- nqi -3x3 matrix in Hz, ready for input
- into create.m function

## Implementation structure

- Converts CASTEP EFG tensor (it is printed in atomic units) to NQI
- 3x3 tensor in Hz that is required by Spinach. Syntax:
- nqi=castep2nqi(V,Q,I)
- V -EFG tensor from CASTEP output, a.u.
- Q -nuclear quadrupole moment, barn
- I -nuclear spin quantum number
- nqi -3x3 matrix in Hz, ready for input
- into create.m function
- Check consistency
- Fundamental constants
- Calculation
- Consistency enforcement
