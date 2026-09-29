# kernel/conventions/transforms/ang2cgsppm.m

[Source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/ang2cgsppm.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=ang2cgsppm.m)

## Contract

ang2cgsppm defines a scalar unit conversion for magnetic susceptibility values used by Spinach pseudocontact-shift functionality: values expressed in cubic Angstrom are converted to cgs-ppm, also quoted as cm^3/mol by quantum-chemistry packages. It does not rotate a tensor or change its array shape.

## Inputs and output

- ang: any numeric array of susceptibility values in cubic Angstrom. The source checks numeric type; it does not impose a tensor shape or a real-only restriction.
- cgsppm: the converted array, with the same element arrangement as the input and units cgs-ppm (cm^3/mol).

The source formula is cgsppm = 6.02214129e23 * ang / (4*pi*1e18). The numeric factor is part of the implementation's stated conversion and is preserved here.

## Source-supported use

The documented call is cgsppm=ang2cgsppm(ang), where ang is an array in cubic Angstrom.
