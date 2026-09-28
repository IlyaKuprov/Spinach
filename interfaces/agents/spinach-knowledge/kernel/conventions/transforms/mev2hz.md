# kernel/conventions/transforms/mev2hz.m

- Signature: `hz=mev2hz(mev)`

## Purpose

Converts energy values from millielectronvolts (meV) to frequency in hertz (Hz).

## Physical / mathematical content

Using `E=h*nu` and the exact electronvolt-to-joule conversion, `hz=(1e-3*1.602176634e-19/6.62607015e-34)*mev`.

## Numerical / algorithmic content

The conversion is a constant elementwise scaling; there are no optional arguments or defaults.

## Syntax

```matlab
hz=mev2hz(mev)
```

## Parameters / inputs

- `mev` — a real numeric array of any dimensions, containing energies in meV.

## Outputs

- `hz` — an array of frequencies in Hz with the same dimensions as `mev`.

## Implementation structure

The function checks that the input is numeric and real, then applies the conversion factor.
