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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 32-33: Check consistency; implemented by `grumble(angles1,weights1,angles2,weights2)`.
- Lines 35-36: Convert both grids to quaternions; implemented by `qter1=euler2qter(angles1(:,1),angles1(:,2),angles1(:,3))`.
- Lines 39-40: Build a table of quaternion pairs; implemented by `n1=size(angles1,1); n2=size(angles2,1)`.
- Lines 46-47: Multiply up quaternions; implemented by `qter.u=u1.*u2-i1.*i2-j1.*j2-k1.*k2`.
- Lines 52-53: Convert quaternions into angles; implemented by `[alphas,betas,gammas]=qter2euler(qter); angles=[alphas betas gammas]`.
- Lines 55-56: Tile the weights; implemented by `weights=kron(weights1,weights2)`.

### Key state/data transformations

- Lines 36: computes `qter1` using `qter1=euler2qter(angles1(:,1),angles1(:,2),angles1(:,3))`.
- Lines 37: computes `qter2` using `qter2=euler2qter(angles2(:,1),angles2(:,2),angles2(:,3))`.
- Lines 40: computes `n1` using `n1=size(angles1,1); n2=size(angles2,1)`.
- Lines 41: computes `u1` using `u1=kron(qter1.u,ones(n2,1)); u2=kron(ones(n1,1),qter2.u)`.
- Lines 42: computes `i1` using `i1=kron(qter1.i,ones(n2,1)); i2=kron(ones(n1,1),qter2.i)`.
- Lines 43: computes `j1` using `j1=kron(qter1.j,ones(n2,1)); j2=kron(ones(n1,1),qter2.j)`.
- Lines 44: computes `k1` using `k1=kron(qter1.k,ones(n2,1)); k2=kron(ones(n1,1),qter2.k)`.
- Lines 47: computes `qter.u` using `qter.u=u1.*u2-i1.*i2-j1.*j2-k1.*k2`.
- Lines 48: computes `qter.i` using `qter.i=u1.*i2+i1.*u2+j1.*k2-k1.*j2`.
- Lines 49: computes `qter.j` using `qter.j=u1.*j2-i1.*k2+j1.*u2+k1.*i2`.
- Lines 50: computes `qter.k` using `qter.k=u1.*k2+i1.*j2-j1.*i2+k1.*u2`.
- Lines 53: computes `[alphas,betas,gammas]` using `[alphas,betas,gammas]=qter2euler(qter); angles=[alphas betas gammas]`.
- Lines 56: computes `weights` using `weights=kron(weights1,weights2)`.

### Local helper functions

- Line 61: `grumble()` — `function grumble(angles1,weights1,angles2,weights2)`.
  - Representative operation: `if (~isnumeric(angles1))||(~isreal(angles1))||(size(angles1,2)~=3)`.
  - Representative operation: `error('angles1 must be a real matrix with three columns.')`.

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
