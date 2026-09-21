# interfaces/pdb_bmrb/read_pdb_nuc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/pdb_bmrb/read_pdb_nuc.m`
- Signature: `[res_num,res_typ,pdb_id,coords]=read_pdb_nuc(pdb_file_name)`
- Total lines: 95

## Purpose

Reads the coordinates of all atoms from the user-specified PDB file and returns, for each atom, the residue number, the residue type, the PDB label and the Cartesian coordinates. Syntax: [res_num,res_typ,pdb_id,coords]=read_pdb_nuc(pdb_file_name)

## Physical / mathematical content

- PDB/BMRB interfaces. These files bridge biomolecular structure/assignment data and Spinach input structures, including atom selection, coordinates, and chemical-shift metadata.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- pdb_file_name -a character string

## Outputs

- nuc_num -nspins x 1 vector giving the number of the
- nucleotide to which each spin belongs
- nuc_typ -nspins x 1 cell array of strings giving the
- PDB identifier of the nucleotide to which
- each spin belongs (e.g. 'GUA')
- pdb_id -nspins x 1 cell array of strings giving the
- PDB identifier of the nucleic acid atom
- type to which each spin belongs (e.g. 'C1P')
- coords -nspins x 1 cell array of 3-vectors giving
- Cartesian coordinates of each spin in Angstrom
- Note: All atoms in the file are read, make sure the PDB only contains
- one model and one chain. Chain identifiers are accepted but not
- returned, files with more than one chain are refused.

## Implementation structure

- Reads the coordinates of all atoms from the user-specified PDB file
- and returns, for each atom, the residue number, the residue type, the
- PDB label and the Cartesian coordinates. Syntax:
- [res_num,res_typ,pdb_id,coords]=read_pdb_nuc(pdb_file_name)
- pdb_file_name -a character string
- nuc_num -nspins x 1 vector giving the number of the
- nucleotide to which each spin belongs
- nuc_typ -nspins x 1 cell array of strings giving the
- PDB identifier of the nucleotide to which
- each spin belongs (e.g. 'GUA')
- pdb_id -nspins x 1 cell array of strings giving the
- PDB identifier of the nucleic acid atom

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fopen()`, `feof()`, `fgetl()`, `textscan()`, `all()`, `cellfun()`, `res_num()`, `upper()`, `fclose()`, `ischar()`.
