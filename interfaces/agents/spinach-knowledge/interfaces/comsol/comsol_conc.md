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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 25-26: Check consistency; implemented by `grumble(spin_system,file_name)`.
- Lines 28-29: Open the file; implemented by `fid=fopen(file_name,'r')`.
- Lines 31-32: Concentration readout count; implemented by `while true`.
- Lines 42-43: Read concentrations; implemented by `X=nan(ncon,1); Y=nan(ncon,1); C=cell(ncon,1)`.
- Lines 51-52: Close the file; implemented by `fclose(fid)`.
- Lines 54-56: Check that vertex locations are the same; implemented by `if (norm(spin_system.mesh.x-X,1)>1e-6)|| (norm(spin_system.mesh.y-Y,1)>1e-6)`.

### Control flow inferred from the code

- Line 32: `while` loop over `true`.
- Line 34: conditional branch on `contains(next_line,'% Nodes:')`.
- Line 44: `for` loop over `n=1:ncon`.
- Line 55: conditional branch on `(norm(spin_system.mesh.x-X,1)>1e-6)||`.

### Key state/data transformations

- Lines 29: computes `fid` using `fid=fopen(file_name,'r')`.
- Lines 33: computes `next_line` using `next_line=fgetl(fid)`.
- Lines 35: computes `ncon` using `ncon=textscan(next_line,'%s %s %f')`.
- Lines 43: computes `X` using `X=nan(ncon,1); Y=nan(ncon,1); C=cell(ncon,1)`.
- Lines 45: computes `VS` using `VS=textscan(fgetl(fid),'%f'); VS=VS{1}`.
- Lines 46: computes `X(n)` using `X(n)=VS(1); Y(n)=VS(2); C{n}=VS(4:end)'`.
- Lines 48: computes `spin_system.mesh.c` using `spin_system.mesh.c=cell2mat(C)`.

### Local helper functions

- Line 63: `grumble()` — `function grumble(spin_system,file_name)`. Feed a captive wolf as much as you like -he would still leg it into the forest at the first opportunity.
  - Representative operation: `if ~ischar(file_name)`.
  - Representative operation: `error('file_name must be a character string.')`.

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
