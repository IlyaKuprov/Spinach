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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(V,Q,I)`.
- Lines 28-29: Fundamental constants; implemented by `efg_atomic=9.717362e+21`.
- Lines 33-34: Calculation; implemented by `nqi=V*efg_atomic*(Q*1e-28)*e_charge/(h_planck*2*I*(2*I-1))`.

### Key state/data transformations

- Lines 29: computes `efg_atomic` using `efg_atomic=9.717362e+21`.
- Lines 30: computes `e_charge` using `e_charge=1.60217657e-19`.
- Lines 31: computes `h_planck` using `h_planck=6.62606957e-34`.
- Lines 34: computes `nqi` using `nqi=V*efg_atomic*(Q*1e-28)*e_charge/(h_planck*2*I*(2*I-1))`.

### Local helper functions

- Line 39: `grumble()` — `function grumble(V,Q,I)`.
  - Representative operation: `if (~isnumeric(V))||(~isnumeric(Q))||(~isnumeric(I))`.
  - Representative operation: `error('all inputs must be numeric.')`.

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
