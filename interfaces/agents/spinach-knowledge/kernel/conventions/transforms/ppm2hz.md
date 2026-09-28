# kernel/conventions/transforms/ppm2hz.m

- Signature: `hz=ppm2hz(ppm,B0,nucleus)`

## Purpose

Converts a chemical shift in parts per million (ppm) to a resonance offset in hertz (Hz).

## Physical / mathematical content

The source computes `hz=1e-6*ppm*(B0*spin(nucleus)/(2*pi))`; the sign of `spin(nucleus)` is retained.

## Numerical / algorithmic content

All three arguments are required. The conversion is linear in `ppm`, with the nucleus-dependent gyromagnetic ratio supplied by `spin(nucleus)`.

## Syntax

```matlab
hz=ppm2hz(ppm,B0,nucleus)
```

## Parameters / inputs

- `ppm` — real numeric chemical-shift value or array in ppm.
- `B0` — real numeric scalar magnetic induction in tesla.
- `nucleus` — character array specifying the isotope, such as `'1H'`.

## Outputs

- `hz` — resonance offset in Hz, with the sign determined by `spin(nucleus)`.

## Implementation structure

The function validates the arguments and uses the specified isotope's spin conversion in the frequency calculation.
