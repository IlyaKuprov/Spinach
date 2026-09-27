# interfaces/comsol/comsol_conc.m

- Signature: `spin_system=comsol_conc(spin_system,file_name)`

## Purpose

Imports ASCII 2D concentration files produced by COMSOL. Syntax: spin_system=comsol_conc(spin_system,file_name)

## Physical / mathematical content

Imports concentration values associated with vertices of an existing Spinach mesh. The file vertex coordinates must match `spin_system.mesh.x` and `spin_system.mesh.y`.

## Numerical / algorithmic content

Finds the line containing `% Nodes:` and reads its third field as the concentration-record count. After skipping four lines, reads that many coordinate and concentration records, assigns the concentration data to `spin_system.mesh.c`, then checks the imported coordinates against the mesh with a 1-norm tolerance of `1e-6`.

## Parameters / inputs

- spin_system -Spinach spin system object
- file_name -a character string

## Outputs

- the following fields are added to spin_system object
- mesh.c -stack of column vectors with
- concentrations (in rows) at
- each vertex of the mesh

## Implementation structure

Validates that `file_name` is a character string and that `spin_system.mesh` exists; opens the file; finds `% Nodes:` and reads the record count; reads the vertex coordinates and concentrations; assigns `spin_system.mesh.c`; closes the file; and checks the vertex coordinates against the existing mesh.
