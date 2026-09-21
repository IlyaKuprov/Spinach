# kernel/grids/grid_kron.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/kernel/grids/grid_kron.m`
- Signature: `[angles,weights]=grid_kron(angles1,weights1,angles2,weights2)`
- Total lines: 90

## Purpose

Spherical grid direct product. Tiles one grid using the rotations of the other. Grids should be supplied using Euler angles in three col- umns [alphas betas gammas] in radians. Syntax: [angles,weights]=grid_kron(angles1,weights1,angles2,weights2)

## Physical / mathematical content

- Quadrature and geometry utilities. These files generate spherical/SO(3) grids, Voronoi weights, and adaptive integration tools for orientation averaging.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.
- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- angles1 -Euler angles (ZYZ active) of the first grid,
- as [alpha beta gamma], radians
- weights1 -weights of the first grid
- angles2 -Euler angles (ZYZ active) of the second grid,
- as [alpha beta gamma], radians
- weights1 -weights of the second grid

## Outputs

- angles -Euler angles (ZYZ active) of the product grid,
- as [alpha beta gamma], rad
- weights -weights of the product grid

## Implementation structure

- Spherical grid direct product. Tiles one grid using the rotations of
- the other. Grids should be supplied using Euler angles in three col-
- umns [alphas betas gammas] in radians. Syntax:
- [angles,weights]=grid_kron(angles1,weights1,angles2,weights2)
- angles1 -Euler angles (ZYZ active) of the first grid,
- as [alpha beta gamma], radians
- weights1 -weights of the first grid
- angles2 -Euler angles (ZYZ active) of the second grid,
- weights1 -weights of the second grid
- angles -Euler angles (ZYZ active) of the product grid,
- as [alpha beta gamma], rad
- weights -weights of the product grid

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `euler2qter()`, `angles1()`, `angles2()`, `qter2euler()`, `any()`.
