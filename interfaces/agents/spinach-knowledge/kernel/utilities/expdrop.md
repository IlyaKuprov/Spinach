# kernel/utilities/expdrop.m

- Signature: `drop=expdrop(from,to,duration,npoints,drop_rate)`

## Purpose

Exponential drop function. Produces an exponential fall-off from a specified value to a specified value with the specified rate and the number of points. Syntax: drop=expdrop(from,to,duration,npoints,drop_rate)

## Physical / mathematical content

- Constructs a sampled exponential decay curve that starts at `from` and reaches `to` at `duration`, using `drop_rate` and `npoints` samples.

## Numerical / algorithmic content

## Parameters / inputs

- from -the value to drop from
- to -the value to drop to
- duration -drop duration, seconds
- npoints -the number of discretisation points
- in the drop
- drop_rate -exponential drop rate, Hz

## Outputs

- drop -a row vector with the fall-off

## Implementation structure

- Exponential drop function. Produces an exponential fall-off
- from a specified value to a specified value with the specified
- rate and the number of points. Syntax:
- drop=expdrop(from,to,duration,npoints,drop_rate)
- from -the value to drop from
- to -the value to drop to
- duration -drop duration, seconds
- npoints -the number of discretisation points
- in the drop
- drop_rate -exponential drop rate, Hz
- drop -a row vector with the fall-off
- Check consistency
