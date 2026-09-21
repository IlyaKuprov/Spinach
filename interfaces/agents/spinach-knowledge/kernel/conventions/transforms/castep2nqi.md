# kernel/conventions/transforms/castep2nqi.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/conventions/transforms/castep2nqi.m`
- Signature: `nqi=castep2nqi(V,Q,I)`
- Total lines: 64

## Purpose

Converts CASTEP EFG tensor (it is printed in atomic units) to NQI 3x3 tensor in Hz that is required by Spinach. Syntax: nqi=castep2nqi(V,Q,I)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `all()`.
