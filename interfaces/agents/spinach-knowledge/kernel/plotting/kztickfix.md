# kernel/plotting/kztickfix.m

- Signature: `kztickfix()`

## Purpose

Switches Z axis tick labels to engineering notation by setting the numeric-ruler exponent to a multiple of three and installs a pan/zoom auto-updater. Syntax: kztickfix()

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- none

## Outputs

- updates the current axis system and installs an auto-updater

## Implementation structure

- Switches Z axis tick labels to engineering notation by setting
- the numeric-ruler exponent to a multiple of three and installs
- a pan/zoom auto-updater. Syntax:
- kztickfix()
- none
- updates the current axis system and installs an auto-updater
- Current Z axis ruler
- Check consistency
- Preserve any existing limit-change callback
- Install the automatic exponent updater
- Set the current exponent
- Local update callback
