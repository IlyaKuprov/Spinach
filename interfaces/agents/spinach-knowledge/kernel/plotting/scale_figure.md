# kernel/plotting/scale_figure.m

- Signature: `scale_figure(by)`

## Purpose

Scales the current figure from the default size by the factors provided by the user. Syntax: scale_figure(by)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- by -two-element row vector of scaling factors,
- format: [width height]

## Implementation structure

- Scales the current figure from the default size by the factors
- provided by the user. Syntax:
- scale_figure(by)
- by - two-element row vector of scaling factors,
- format: [width height]
- Check consistency
- Get figure location
- Get figure centroid
- Get default figure size
- Scale the figure
- Update figure parameters
- Consistency enforcement
