# kernel/grids/arclength.m

- Signature: `sig=arclength(r1,r2)`

## Purpose

Arc length between two points on the unit sphere specified by the unit vectors supplied. Syntax: sig=arclength(r1,r2)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

## Parameters / inputs

- r1,r2 -three-element unit vectors with Cartesian
- coordinates of the arc endpoints

## Outputs

- sig -arc length

## Implementation structure

- Arc length between two points on the unit sphere specified by
- the unit vectors supplied. Syntax:
- sig=arclength(r1,r2)
- r1,r2 -three-element unit vectors with Cartesian
- coordinates of the arc endpoints
- sig -arc length
- Check consistency
- Normalise the vectors
- Get the arc length
- Consistency enforcement
- Ah, there's nothing more exciting than science. You get all the
- fun of sitting still, being quiet, writing down numbers, paying
