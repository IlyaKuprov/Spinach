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

- Switches Y axis tick labels to engineering notation by setting
- the numeric-ruler exponent to a multiple of three and installs
- a pan/zoom auto-updater. Syntax:
- kytickfix()
- none
- updates the current axis system and installs an auto-updater
- Current Y axis ruler
- Check consistency
- Preserve any existing limit-change callback
- Install the automatic exponent updater
- Set the current exponent
- Local update callback
