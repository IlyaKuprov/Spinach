# kernel/conventions/transforms/hz2icm.m

- Signature: `icm=hz2icm(hz)`

## Purpose

Converts Hz units used in magnetic resonance into cm^-1 units used in spectroscopy. Syntax: icm=hz2icm(hz) Arrays of any dimensions are supported. Parameters: hz -an array of values in Hz

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Outputs

- icm -an array of values in inverse centimetres

## Implementation structure

- Converts Hz units used in magnetic resonance into cm^-1 units
- used in spectroscopy. Syntax:
- icm=hz2icm(hz)
- Arrays of any dimensions are supported. Parameters:
- hz -an array of values in Hz
- icm -an array of values in inverse centimetres
- Check consistency
- Run the conversion
- Consistency enforcement
- "Damn, would I have to be nice to everyone for two years?!"
- IK, upon being informed that he was to
- organise the 48th ESR Group Conference
