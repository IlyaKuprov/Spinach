# kernel/grids/one_vcell_solidangle.m

- Signature: `S=one_vcell_solidangle(v,centre)`

## Purpose

Solid angle of a convex spherical polygon as described in

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

## Syntax

```matlab
A=one_vcell_solidangle(v,centre)
```

## Parameters / inputs

- v -(3 x n) matrix of unit vectors giving
- the coordinates of each vertex
- centre -centre vertex coordinates, optional

## Outputs

- S -the solid angle, radians

## Implementation structure

- Solid angle of a convex spherical polygon as described in
- A=one_vcell_solidangle(v,centre)
- v -(3 x n) matrix of unit vectors giving
- the coordinates of each vertex
- centre -centre vertex coordinates, optional
- S -the solid angle, radians
- Check consistency
- Straightforward math
- Consistency enforcement
- Being a mathematician is a bit like being a manic
- depressive: you spend your life alternating between
- giddy elation and black despair.
