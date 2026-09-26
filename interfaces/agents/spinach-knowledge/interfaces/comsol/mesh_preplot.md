# interfaces/comsol/mesh_preplot.m

- Signature: `mesh=mesh_preplot(mesh)`

## Purpose

Mesh preprocessing for drawing. Creates edge, triangle, and rectangle data structures needed for fast plotting later.

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Syntax

```matlab
mesh=mesh_preplot(mesh)
```

## Parameters / inputs

- mesh -Spinach mesh object

## Outputs

- mesh -updated mesh object

## Implementation structure

- Mesh preprocessing for drawing. Creates edge, triangle, and
- rectangle data structures needed for fast plotting later.
- mesh=mesh_preplot(mesh)
- mesh -Spinach mesh object
- mesh -updated mesh object
- Check consistency
- Prepare edge array for plotting
- Prepare triangle array for plotting
- Prepare rectangle array for plotting
- Prepare Voronoi cell array for plotting
- Consistency enforcement
- Life is a tragedy for those who feel,
