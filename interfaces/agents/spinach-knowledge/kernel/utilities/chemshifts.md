# kernel/utilities/chemshifts.m

- Signature: `[cs_ppm,cs_hz]=chemshifts(spin_system)`

## Purpose

Returns the chemical shifts of every spin relative to the carrier frequency in the current magnet, in parts per million and hertz.

## Physical / mathematical content

For each spin, the function obtains the isotropic Zeeman frequency from one third of the trace of its Zeeman matrix, then subtracts the corresponding carrier frequency in `spin_system.inter.basefrqs`.

## Numerical / algorithmic content

For each spin `n`, the carrier-subtracted value `iso` is converted to ppm as `1e6*iso/basefrqs(n)` and to hertz as `-iso/(2*pi)`.

## Parameters / inputs

- `spin_system` — spin system descriptor object

## Outputs

- `cs_ppm` — chemical shifts in ppm
- `cs_hz` — chemical shifts in Hz
