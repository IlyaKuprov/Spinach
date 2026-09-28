# kernel/plotting/kytickfix.m

- Signature: `kytickfix()`

## Purpose

Switches Y axis tick labels to engineering notation by setting the numeric-ruler exponent to a multiple of three and installs a pan/zoom auto-updater. Syntax: kytickfix()

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- none

## Outputs

- updates the current axis system and installs an auto-updater

## Implementation structure

- Uses the current axes Y ruler and rejects non-numeric or non-linear rulers.
- Installs a limits-changed callback while preserving and later invoking any existing callback.
- Sets the ruler exponent from finite, non-zero tick values (or limits when no tick values are available) to a multiple of three; if neither provides values, it sets the exponent to zero.

[Source reference](https://spindynamics.org/wiki/index.php?title=kytickfix.m)
