# kernel/kinetics/flow_gen.m

- Signature: `F=flow_gen(spin_system,parameters)`

## Purpose

Builds a spatial motion generator for hydrodynamic flow and diffusion on the Voronoi mesh stored in `spin_system.mesh`.

## Physical / mathematical content

The generator represents transfers between adjacent Voronoi cells. For each shared cell boundary, the code estimates the advective contribution from the boundary length, the distance between cell centres, and the average velocity in the two cells. Diffusion contributes transfers proportional to `parameters.diff`; diagonal entries are set to balance the off-diagonal transfers, and Voronoi cell areas are applied to the generator.

## Numerical / algorithmic content

The routine identifies neighboring cells through mesh triangles and shared Voronoi vertices, assembles local transfer entries in a `parfor` loop, then constructs a sparse matrix. If `parameters.diff` is absent, it defaults to zero.

## Parameters / inputs

- `spin_system` - Spinach system descriptor with mesh data, including Voronoi cells, active mesh vertices, coordinates, and velocity components; mesh subfields are produced by COMSOL import functions.
- `parameters.diff` - diffusion coefficient in m^2/s (default: zero).

## Outputs

- `F` - spatial motion generator matrix with dimension equal to the number of Voronoi cells.

## Implementation structure

The routine validates the mesh and its indexing and Voronoi data, finds neighboring cells sharing a Voronoi edge, computes advection and diffusion transfers, and assembles and balances the sparse generator.
