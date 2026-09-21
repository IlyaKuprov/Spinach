# interfaces/comsol/comsol_conc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/comsol/comsol_conc.m`
- Signature: `spin_system=comsol_conc(spin_system,file_name)`
- Total lines: 76

## Purpose

Imports ASCII 2D concentration files produced by COMSOL. Syntax: spin_system=comsol_conc(spin_system,file_name)

## Physical / mathematical content

- COMSOL interfaces. These files are mostly data-structure and numerical-geometry utilities for bringing concentration, velocity, and mesh data from finite-element simulations into Spinach transport calculations.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- spin_system -Spinach spin system object
- file_name -a character string

## Outputs

- the following fields are added to spin_system object
- mesh.c -stack of column vectors with
- concentrations (in rows) at
- each vertex of the mesh

## Implementation structure

- Imports ASCII 2D concentration files produced by COMSOL. Syntax:
- spin_system=comsol_conc(spin_system,file_name)
- spin_system -Spinach spin system object
- file_name -a character string
- the following fields are added to spin_system object
- mesh.c -stack of column vectors with
- concentrations (in rows) at
- each vertex of the mesh
- Check consistency
- Open the file
- Concentration readout count
- Read concentrations

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fopen()`, `fgetl()`, `contains()`, `textscan()`, `report()`, `num2str()`, `nan()`, `cell2mat()`, `fclose()`, `ischar()`, `isfield()`.
