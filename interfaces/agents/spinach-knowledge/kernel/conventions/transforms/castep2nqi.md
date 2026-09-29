# kernel/conventions/transforms/castep2nqi.m

- Signature: `nqi=castep2nqi(V,Q,I)`

## Purpose

Converts a CASTEP electric-field-gradient tensor into the nuclear quadrupole interaction tensor used by Spinach.

## Conversion

For a real 3-by-3 EFG tensor `V` in atomic units, nuclear quadrupole moment `Q` in barns, and spin quantum number `I`, the implemented conversion is

`nqi = V * 9.717362e21 * (Q * 1e-28) * 1.60217657e-19 / (6.62606957e-34 * 2 * I * (2 * I - 1))`.

The scalar factor combines the source's EFG atomic-unit constant, barn-to-square-metre factor, elementary charge, Planck constant, and spin denominator. It scales `V` without a coordinate rotation or tensor reparameterisation, so the output is a 3-by-3 tensor in Hz, ready for `create.m`.

## Inputs and outputs

- `V`: real numeric 3-by-3 CASTEP EFG tensor in atomic units.
- `Q`: real numeric scalar nuclear quadrupole moment in barns.
- `I`: real numeric scalar integer or half-integer, at least 1.
- `nqi`: 3-by-3 nuclear quadrupole interaction tensor in Hz.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/castep2nqi.m)
- [Spinach Wiki: castep2nqi.m](https://spindynamics.org/wiki/index.php?title=castep2nqi.m)
