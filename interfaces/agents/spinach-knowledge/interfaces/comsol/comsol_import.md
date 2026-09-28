# interfaces/comsol/comsol_import.m

- Signature: `mesh=comsol_import(comsol)`

## Purpose

COMSOL 2D mesh data import, cropping and preprocessing for Spinach. Syntax: mesh=comsol_import(comsol)

## Physical / mathematical content

Imports a two-dimensional COMSOL mesh and vertex-centred flow velocities into a Spinach mesh object, then crops the mesh, inactivates specified vertices, computes a Voronoi tessellation, and prepares plotting data.

## Numerical / algorithmic content

Calls `comsol_mesh`, `comsol_velo`, `mesh_crop`, `mesh_inact`, `mesh_vorn`, and `mesh_preplot` in that order.

## Parameters / inputs

- comsol.mesh_file -name of an ASCII file with
- vertex coordinates and edge
- index produced by COMSOL
- comsol.velo_file -name of an ASCII file with
- vertex-centred flow veloci-
- ties produced by COMSOL
- comsol.crop -{[xmin xmax],[ymin ymax]}
- region of the mesh to retain
- comsol.inactivate -a row vector with mesh vertex
- indices to deactivate

## Outputs

- mesh – Spinach mesh object with
- ▸ geometry (.x, .y, .idx)
- ▸ flow velocities (.u, .v)
- ▸ Voronoi tessellation (.vor)
- ▸ fast-plot auxiliaries (.plot)
- Notes: internally this routine is just a convenience wrapper
- that calls the following functions
- ▸ comsol_mesh() – mesh import
- ▸ comsol_velo() – velocity import
- ▸ mesh_crop() – region trimming
- ▸ mesh_vorn() – Voronoi tessellation
- ▸ mesh_preplot() – plotting accelerators

## Implementation structure

Validates the COMSOL input structure and required fields, then imports the mesh and velocities, crops the region, inactivates specified vertices, computes the Voronoi tessellation, and prepares plotting data.
