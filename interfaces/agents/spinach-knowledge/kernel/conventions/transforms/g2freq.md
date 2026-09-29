# kernel/conventions/transforms/g2freq.m

## Purpose

Converts g-values and a scalar magnetic field to electron Zeeman frequency values.

## Signature

`f=g2freq(g,B)`

## Conversion

The source first evaluates `omega = B*spin('E')/(2*pi)`, then `f = g .* omega / 2.0023193043622`. Equivalently:

`f = g .* (B*spin('E')/(2*pi)) / 2.0023193043622`

The source documents `B` in tesla and `f` in Hz. The denominator is the source's free-electron g-factor constant, `2.0023193043622`.

## Inputs and output

- `g`: finite, real numeric array; no positivity or orientation constraint is applied.
- `B`: finite real numeric scalar in tesla.
- `f`: frequency array in Hz, with the same dimensions as `g`.

Only a scalar field is accepted; the function does not take an orientation or apply a tensor rotation. The frequency scale uses `spin('E')` as shown in the source.

## References

- MATLAB source: [kernel/conventions/transforms/g2freq.m](https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/g2freq.m)
- [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=g2freq.m)
