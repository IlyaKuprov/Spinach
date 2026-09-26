# examples/microfluidics/show_mesh.m

- Signature: `show_mesh()`

## Purpose

Import, Voronoi tessellation, and plotting of the hydrodynamic mesh and velocity field from COMSOL.

## Physical / mathematical content

- Microfluidics examples. The coupled model is spin dynamics plus advection-diffusion-reaction transport on a mesh or regular grid. Numerical issues include finite-difference operators, mesh interpolation, and coupled reaction-flow evolution.

## Numerical / algorithmic content

- Numerical integration over angles or geometry is part of the implementation, so point placement and weights are as important as the local Hamiltonian calculations.

## Implementation structure

- Import, Voronoi tessellation, and plotting of the
- hydrodynamic mesh and velocity field from COMSOL.
- Import hydrodynamics information
- No spin system here
- Draw the mesh
