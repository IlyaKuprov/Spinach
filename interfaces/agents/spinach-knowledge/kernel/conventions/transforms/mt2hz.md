# kernel/conventions/transforms/mt2hz.m

Source implementation: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/mt2hz.m
Spinach Wiki: [mt2hz.m](https://spindynamics.org/wiki/index.php?title=mt2hz.m)

## Purpose

Converts hyperfine couplings from milliTesla (mT) to hertz (Hz, linear frequency). The field specification is described as the magnetic field at which the electron frequency equals the supplied frequency.

## Syntax

```matlab
hfc_hz=mt2hz(hfc_mt,g)
```

The second argument may be omitted; in that case the function displays a message and uses the free-electron g-factor `2.0023193043622`.

## Inputs

- `hfc_mt` is a numeric array of real values in mT; arrays of any dimensions are accepted.
- `g` is a real numeric scalar. If provided for an isotope-specific calculation, use the desired g-factor explicitly; only the omitted-argument case receives the free-electron default.

## Output and conversion

- `hfc_hz` is an array of values in Hz with the same dimensions as `hfc_mt`.

The implementation defines `muB=9.274009994*10^-24`, `hbar=1.054571628*10^-34`, and `C=1e-3*g*muB/(hbar*2*pi)`, then returns `hfc_hz=C*hfc_mt`. These source constants and the mT factor `1e-3` are retained as implemented.

Validation checks that `hfc_mt` is numeric and real, and that `g` is numeric, real, and scalar. The implementation does not impose an additional positivity or finiteness check on `g`.
