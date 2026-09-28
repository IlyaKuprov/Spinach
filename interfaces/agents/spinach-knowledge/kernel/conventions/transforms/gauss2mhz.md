# kernel/conventions/transforms/gauss2mhz.m

- Signature: `hfc_mhz=gauss2mhz(hfc_gauss,g)`

## Purpose

Converts hyperfine couplings from gauss to MHz (linear frequency). Here, a gauss value can be specified as the magnetic field at which the electron frequency equals the frequency provided.

## Parameters / inputs

- `hfc_gauss`: real numeric array of hyperfine couplings in gauss; arrays of any dimensions are supported.
- `g`: optional real scalar electron g-factor. If omitted, the free-electron value `2.0023193043622` is used.

## Output

- `hfc_mhz`: array of values in MHz, with the same shape as `hfc_gauss`.

## Conversion

The conversion uses `hfc_mhz = 1e-10 * g * muB * hfc_gauss / (hbar * 2*pi)`, with the Bohr magneton and reduced Planck constant in SI units.

Source: [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=gauss2mhz.m)
