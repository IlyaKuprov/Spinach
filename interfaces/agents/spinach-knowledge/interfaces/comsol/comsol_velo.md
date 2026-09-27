# interfaces/comsol/comsol_velo.m

- Signature: `mesh=comsol_velo(mesh,file_name)`

## Purpose

Imports ASCII 2D flow velocity files produced by COMSOL. Syntax: mesh=comsol_velo(mesh,file_name)

## Physical / mathematical content

The imported `mesh.u` and `mesh.v` values are velocity components at mesh vertices. The velocity-file vertex coordinates must match `mesh.x` and `mesh.y`.

## Numerical / algorithmic content

Requires `file_name` to be a character string. Reads the velocity-readout count from the `% Nodes:` line, then parses each five-column row as vertex coordinates and velocity components. Rejects coordinates that differ from `mesh.x` or `mesh.y` by more than `1e-6` in 1-norm; otherwise stores the velocities in `mesh.u` and `mesh.v`.

## Parameters / inputs

- mesh -mesh object produced by comsol_mesh()
- file_name -a character string

## Outputs

- the following fields are added to the mesh object
- mesh.u, mesh.v -column vectors with velocities
- at each vertex of the mesh

## Implementation structure

Validates `file_name`, opens the file, finds `% Nodes:` to determine the readout count, parses coordinates and velocity components, and closes the file. It checks the coordinates against the mesh before assigning `mesh.u` and `mesh.v`.
