# kernel/conventions/transforms/castep2nqi.m

- Signature: `nqi=castep2nqi(V,Q,I)`

## Purpose

Converts CASTEP EFG tensor (it is printed in atomic units) to NQI 3x3 tensor in Hz that is required by Spinach. Syntax: nqi=castep2nqi(V,Q,I)

## Physical / mathematical content
The conversion multiplies the CASTEP EFG tensor V by a scalar factor to obtain a 3×3 nuclear quadrupole interaction tensor in Hz. The factor converts EFG atomic units and Q in barns to SI units, then applies the nuclear charge and divides by Planck's constant and 2I(2I−1). No coordinate rotation or tensor reparameterisation is performed.

## Numerical / algorithmic content
The calculation is `nqi=V*9.717362e21*(Q*1e-28)*1.60217657e-19/(6.62606957e-34*2*I*(2*I-1))`. The result retains V's 3×3 shape. Before calculation, the function requires real numeric inputs, a 3×3 V, a scalar Q, and a scalar integer or half-integer I of at least 1.

## Parameters / inputs

- V -EFG tensor from CASTEP output, a.u.
- Q -nuclear quadrupole moment, barn
- I -nuclear spin quantum number

## Outputs

- nqi -3x3 matrix in Hz, ready for input
- into create.m function

## Implementation structure
The main function calls the local `grumble(V,Q,I)` validator, defines the EFG atomic-unit conversion factor, elementary charge, and Planck constant, then computes `nqi` by scalar multiplication of V. `grumble` raises errors for nonnumeric or nonreal inputs, an incorrectly sized V, a nonscalar Q, or an I that is not a scalar integer or half-integer of at least 1.
