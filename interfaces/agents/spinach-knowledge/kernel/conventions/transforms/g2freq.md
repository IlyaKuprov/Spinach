# kernel/conventions/transforms/g2freq.m

- Signature: `f=g2freq(g,B)`

## Purpose

Converts electron `g`-values and magnetic field to electron Zeeman frequencies.

## Parameters / inputs

- `g`: finite real scalar or numeric array of g-values.
- `B`: finite real scalar magnetic field in tesla.

## Output

- `f`: frequency in Hz, with the same array shape as `g`.

## Conversion

`f = g * B * spin('E') / (2*pi*2.0023193043622)`

The denominator uses the free-electron g-factor `2.0023193043622` to scale the free-electron carrier frequency returned by `spin('E')`.

Source: [Spinach Wiki](https://spindynamics.org/wiki/index.php?title=g2freq.m)
