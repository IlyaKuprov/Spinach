# kernel/plotting/bwr_cmap.m

- Signature: `cmap=bwr_cmap()`

## Purpose

Builds a 255-by-3 RGB blue–white–red colour map, with white at the midpoint for zero-valued data.

## Numerical / algorithmic content

The first 128 rows interpolate from blue to white; rows 128–255 interpolate from white to red. The component values are then squared (`cmap=cmap.^2`) to apply the quadratic contrast curve.

## Outputs

- `cmap` — 255-by-3 MATLAB RGB colour map.
