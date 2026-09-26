# kernel/conventions/transforms/mev2hz.m

- Signature: `hz=mev2hz(mev)`

## Purpose

Converts meV energy units used in solid state physics and phonon spectroscopy into Hz units preferred in magnetic resonance. Syntax: hz=mev2hz(mev) Arrays of any dimensions are supported. Parameters: mev -an array of values in milli-electronvolts

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Outputs

- hz -an array of values in Hz

## Implementation structure

- Converts meV energy units used in solid state physics and
- phonon spectroscopy into Hz units preferred in magnetic
- resonance. Syntax:
- hz=mev2hz(mev)
- Arrays of any dimensions are supported. Parameters:
- mev -an array of values in milli-electronvolts
- hz -an array of values in Hz
- Check consistency
- Run the conversion
- Consistency enforcement
- O God, I could be bounded in a nutshell, and count
- myself a king of infinite space, were it not that I
