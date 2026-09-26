# kernel/conventions/transforms/icm2hz.m

- Signature: `hz=icm2hz(icm)`

## Purpose

Converts cm^-1 units used in spectroscopy into Hz units preferred in magnetic resonance. Syntax: hz=icm2hz(icm) Arrays of any dimensions are supported. Parameters: icm -an array of values in inverse centimetres

## Physical / mathematical content

- Convention and tensor-transform utilities. They convert among tensor parameterisations, coordinate systems, and unit systems; the underlying mathematics is linear algebra on rank-2 tensors and rotation representations.

## Numerical / algorithmic content

## Outputs

- hz -an array of values in Hz

## Implementation structure

- Converts cm^-1 units used in spectroscopy into Hz units
- preferred in magnetic resonance. Syntax:
- hz=icm2hz(icm)
- Arrays of any dimensions are supported. Parameters:
- icm -an array of values in inverse centimetres
- hz -an array of values in Hz
- Check consistency
- Run the conversion
- Consistency enforcement
- Satire thrives where the usual checks on human
- folly fail.
- The Economist
