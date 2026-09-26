# kernel/plotting/cylgrid.m

- Signature: `cylgrid(zmin,zmax,rmax)`

## Purpose

Draws a cylindrical grid with 10% spacing added around the indicated data extent values. Syntax: cylgrid(zmin,zmax,rmax)

## Physical / mathematical content

## Numerical / algorithmic content

## Parameters / inputs

- zmin -lower bound on the Z axis
- zmax -upper bound on the Z axis
- rmax -upper bound on the radius

## Outputs

- this function updates the current figure

## Implementation structure

- Draws a cylindrical grid with 10% spacing added around
- the indicated data extent values. Syntax:
- cylgrid(zmin,zmax,rmax)
- zmin -lower bound on the Z axis
- zmax -upper bound on the Z axis
- rmax -upper bound on the radius
- this function updates the current figure
- Check consistency
- Get the extent gaps
- Draw the spokes
- Draw the circles
- Set axis extents
