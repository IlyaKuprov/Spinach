# interfaces/comsol/comsol_velo.m

- Signature: `mesh=comsol_velo(mesh,file_name)`

## Purpose

Imports ASCII 2D flow velocity files produced by COMSOL. Syntax: mesh=comsol_velo(mesh,file_name)

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

## Parameters / inputs

- mesh -mesh object produced by comsol_mesh()
- file_name -a character string

## Outputs

- the following fields are added to the mesh object
- mesh.u, mesh.v -column vectors with velocities
- at each vertex of the mesh

## Implementation structure

- Imports ASCII 2D flow velocity files produced by COMSOL. Syntax:
- mesh=comsol_velo(mesh,file_name)
- mesh -mesh object produced by comsol_mesh()
- file_name -a character string
- the following fields are added to the mesh object
- mesh.u, mesh.v -column vectors with velocities
- at each vertex of the mesh
- Check consistency
- Open the file
- Velocity readout count
- Parse velocity readouts
- Close the file
