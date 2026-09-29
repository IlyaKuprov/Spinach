# kernel/conventions/transforms/hz2ppm.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/hz2ppm.m) · [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=hz2ppm.m)

## Signature

`ppm = hz2ppm(hz, B0, nucleus)`

## Purpose and conversion

Converts a resonance offset to a chemical-shift value with the implemented expression:

`ppm = 1e6 * (2*pi*hz) / (B0*spin(nucleus))`

Here `hz` is the resonance offset in Hz, `B0` is magnetic induction in tesla, and `nucleus` names the isotope (for example, `'1H'`). The source notes that signs of the magnetogyric ratios are preserved. The expression uses `spin(nucleus)` as returned, with no absolute-value operation.

## Inputs and output

- `hz`: real numeric array; no explicit dimensionality or size restriction is checked.
- `B0`: real numeric scalar.
- `nucleus`: character array naming the isotope; the implementation rejects non-character input.
- `ppm`: chemical shift in ppm, returned from the displayed expression.

The source checks that `hz` is numeric and real, that `B0` is numeric, real, and scalar, and that `nucleus` is a character array. It does not explicitly test the sign or nonzero value of `B0`; no explicit size or dimensionality check is applied to `hz`.

## Reference

- [Spinach Wiki: hz2ppm.m](https://spindynamics.org/wiki/index.php?title=hz2ppm.m)
