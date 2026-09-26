# interfaces/comsol/comsol_mesh.m

- Signature: `mesh=comsol_mesh(file_name)`

## Purpose

Imports ASCII 2D mesh files produced by COMSOL. Syntax: mesh=comsol_mesh(file_name)

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

## Parameters / inputs

- file_name -a character string

## Outputs

- mesh.x, mesh.y -column vectors with vertex
- coordinates
- mesh.idx.edges -two-column array of integers
- containing edge index
- mesh.idx.triangles -three-column array of integers
- containing triangle index
- mesh.idx.rectangles -four-column array of integers
- containing rectangle index

## Implementation structure

- Imports ASCII 2D mesh files produced by COMSOL. Syntax:
- mesh=comsol_mesh(file_name)
- file_name -a character string
- mesh.x, mesh.y -column vectors with vertex
- coordinates
- mesh.idx.edges -two-column array of integers
- containing edge index
- mesh.idx.triangles -three-column array of integers
- containing triangle index
- mesh.idx.rectangles -four-column array of integers
- containing rectangle index
- Check consistency
