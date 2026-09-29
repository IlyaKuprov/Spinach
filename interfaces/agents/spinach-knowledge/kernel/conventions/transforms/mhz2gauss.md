# kernel/conventions/transforms/mhz2gauss.m

Source implementation: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/mhz2gauss.m
Spinach Wiki: [mhz2gauss.m](https://spindynamics.org/wiki/index.php?title=mhz2gauss.m)

## Purpose

Converts hyperfine couplings from MHz (linear frequency) to Gauss. The field specification is described as the magnetic field at which the electron frequency equals the supplied frequency.

## Syntax

```matlab
hfc_gauss=mhz2gauss(hfc_mhz,g)
```

The second argument may be omitted; in that case the function displays a message and uses the free-electron g-factor `2.0023193043622`.

## Inputs

- `hfc_mhz` is a numeric array of real values in MHz; arrays of any dimensions are accepted.
- `g` is a real numeric scalar. If provided for an isotope-specific calculation, use the desired g-factor explicitly; only the omitted-argument case receives the free-electron default.

## Output and conversion

- `hfc_gauss` is an array of values in Gauss with the same dimensions as `hfc_mhz`.

The implementation defines `muB=9.274009994*10^-24`, `hbar=1.054571628*10^-34`, and `C=1e-10*g*muB/(hbar*2*pi)`, then returns `hfc_gauss=hfc_mhz/C`. These source constants and the factor `1e-10` are retained as implemented.

Validation checks that `hfc_mhz` is numeric and real, and that `g` is numeric, real, and scalar. The implementation does not impose an additional positivity or finiteness check on `g`.
