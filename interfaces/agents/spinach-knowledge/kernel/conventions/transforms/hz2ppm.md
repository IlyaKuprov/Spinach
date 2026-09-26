# kernel/conventions/transforms/hz2ppm.m

- Signature: `ppm=hz2ppm(hz,B0,nucleus)`

## Purpose

Converts resonance offsets into chemical shifts. Syntax: ppm=hz2ppm(hz,B0,nucleus)

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Parameters / inputs

- hz -resonance offset in Hz
- B0 -magnet induction, Tesla
- nucleus -a string specifying the isotope, e.g. '1H'

## Outputs

- ppm -chemical shift in ppm
- Note: signs of the magnetogyric ratios are preserved.

## Implementation structure

- Converts resonance offsets into chemical shifts. Syntax:
- ppm=hz2ppm(hz,B0,nucleus)
- hz -resonance offset in Hz
- B0 -magnet induction, Tesla
- nucleus -a string specifying the isotope, e.g. '1H'
- ppm -chemical shift in ppm
- Note: signs of the magnetogyric ratios are preserved.
- Check consistency
- Calculate chemical shift in Hz
- Consistency enforcement
- Somebody Else's Wife: oh, I am so tired of being boringly taken
- for granted, he never really acknowledges
