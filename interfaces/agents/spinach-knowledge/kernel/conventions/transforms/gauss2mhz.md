# kernel/conventions/transforms/gauss2mhz.m

## Purpose

Converts hyperfine coupling values in gauss to linear frequency values in MHz. The source describes the Gauss specification as the magnetic field at which the electron frequency equals the supplied frequency.

## Signature

`hfc_mhz=gauss2mhz(hfc_gauss,g)`

## Conversion

The source sets `muB = 9.274009994e-24`, `hbar = 1.054571628e-34`, and computes:

`C = 1e-10 * g * muB / (hbar * 2*pi)`

`hfc_mhz = C * hfc_gauss`

If `g` is omitted, it uses `g = 2.0023193043622` (free-electron g-factor) and displays a message.

## Inputs and output

- `hfc_gauss`: real numeric array in gauss; arrays of any dimensions are supported.
- `g`: optional real numeric scalar. The implementation requires one element but does not explicitly require it to be finite or positive.
- `hfc_mhz`: array in MHz, with the same dimensions as `hfc_gauss`.

The conversion applies a scalar factor; there is no orientation or rotation input.

## References

- MATLAB source: [kernel/conventions/transforms/gauss2mhz.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/gauss2mhz.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=gauss2mhz.m)
