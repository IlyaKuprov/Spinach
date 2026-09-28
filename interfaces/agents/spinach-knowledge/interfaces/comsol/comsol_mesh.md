# interfaces/comsol/comsol_mesh.m

- Signature: `mesh=comsol_mesh(file_name)`

## Purpose

Imports ASCII 2D mesh files produced by COMSOL. Syntax: mesh=comsol_mesh(file_name)

## Physical / mathematical content

Represents a two-dimensional mesh using vertex coordinates and connectivity arrays for edges, triangles, and quadrilaterals.

## Numerical / algorithmic content

Reads the reported mesh-point and element counts from the ASCII file, then parses vertex coordinates and the `edg`, `tri`, and `quad` connectivity sections. Adds one to each element vertex index to convert the file indices to MATLAB indices.

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

Checks that `file_name` is a character string, opens the file, locates the coordinate, `edg`, `tri`, and `quad` sections, reads their counts and data, and closes the file.
