# kernel/conventions/transforms/cgsppm2ang.m

- Signature: `ang=cgsppm2ang(cgsppm)`

## Purpose

Converts magnetic susceptibility values in cgs-ppm (described in the source as cm³/mol) to cubic angstroms for Spinach pseudocontact-shift calculations.

## Conversion

The implementation applies the scalar factor `4*pi*1e18/6.02214129e23`:

`ang = 4*pi*1e18*cgsppm/6.02214129e23`.

The documented output has the same array size as the input and contains susceptibility values in cubic angstroms. The conversion uses the source's numerical Avogadro constant, `6.02214129e23`.

## Inputs and outputs

- `cgsppm`: numeric array of susceptibility values in cgs-ppm (cm³/mol as described in the source). The validator checks numeric type only; it does not impose a real-valued or finite-value condition or a particular shape.
- `ang`: converted susceptibility array, documented as the same size as `cgsppm`, in cubic angstroms.

## References

- [MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/cgsppm2ang.m)
- [Spinach Wiki: cgsppm2ang.m](https://spindynamics.org/wiki/index.php?title=cgsppm2ang.m)
