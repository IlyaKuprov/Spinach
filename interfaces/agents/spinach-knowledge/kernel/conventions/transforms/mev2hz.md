# kernel/conventions/transforms/mev2hz.m

Source implementation: https://github.com/IlyaKuprov/Spinach/blob/main/kernel/conventions/transforms/mev2hz.m
Spinach Wiki: [mev2hz.m](https://spindynamics.org/wiki/index.php?title=mev2hz.m)

## Purpose

Converts energy values in milli-electronvolts (meV), as used in solid-state physics and phonon spectroscopy, to frequency in hertz (Hz), used in magnetic resonance. It applies E = h times nu using the exact electronvolt-to-joule factor in the implementation.

## Syntax

```matlab
hz=mev2hz(mev)
```

## Input

- `mev` is a numeric array of real values in meV. Arrays of any dimensions are accepted; no separate shape restriction is imposed.

## Output and conversion

- `hz` is an array of values in Hz with the same dimensions as `mev`.

The function applies the elementwise scaling

`hz = (1e-3 * 1.602176634e-19 / 6.62607015e-34) * mev`.

The factors are the meV-to-eV multiplier, the exact electronvolt value in joules, and the exact Planck constant in joule-seconds, respectively. The function has no optional arguments or defaults.
