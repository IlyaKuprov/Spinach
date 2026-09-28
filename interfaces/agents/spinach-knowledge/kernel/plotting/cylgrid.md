# kernel/plotting/cylgrid.m

- Signature: `cylgrid(zmin,zmax,rmax)`

## Purpose

Draws a labelled cylindrical grid around the supplied data extent, with a 10% margin in radius and along the z range.

## Parameters / inputs

- `zmin` — lower bound on the z axis; finite real scalar and less than `zmax`.
- `zmax` — upper bound on the z axis; finite real scalar and greater than `zmin`.
- `rmax` — upper bound on the radius; positive finite real scalar.

## Numerical / algorithmic content

The grid adds radial and axial gaps equal to 10% of `rmax` and `zmax-zmin`, respectively. It draws spokes every 30 degrees, concentric circles, and z-axis tick marks and labels, then sets the current axes to a perspective view without default ticks.

## Outputs

Updates the current figure; the function does not return MATLAB outputs.
