# kernel/conventions/transforms/hz2icm.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/hz2icm.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=hz2icm.m)

## Signature

`icm = hz2icm(hz)`

## Purpose and conversion

Converts frequency values used in magnetic resonance to spectroscopic wavenumbers using the implemented equation:

`icm = hz / (100 * 299792458)`

Input `hz` is in Hz; output `icm` is in inverse centimetres (cm^-1). The denominator uses 100 centimetres per metre and the speed-of-light value 299792458 used in the source. Division by this positive constant preserves the sign of the input.

## Inputs and output

- `hz`: real numeric array; arrays of any dimensionality are documented as supported.
- `icm`: converted array in cm^-1, retaining the input array dimensions.

The implementation rejects nonnumeric or nonreal input and imposes no explicit size check.

## Reference

- [Spinach Wiki: hz2icm.m](https://spindynamics.org/wiki/index.php?title=hz2icm.m)
