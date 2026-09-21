# interfaces/pdb_bmrb/read_bmrb.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/pdb_bmrb/read_bmrb.m`
- Signature: `[aa_num,aa_typ,pdb_id,chemsh]=read_bmrb(bmrb_file_name)`
- Total lines: 91

## Purpose

Reads BMRB file, extracts amino acid numbers, amino acid types, PDB atom identifiers and chemical shifts. Syntax: [aa_num,aa_typ,pdb_id,chemsh]=read_bmrb(bmrb_file_name)

## Physical / mathematical content

- PDB/BMRB interfaces. These files bridge biomolecular structure/assignment data and Spinach input structures, including atom selection, coordinates, and chemical-shift metadata.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- bmrb_file_name -a character string with the file name

## Outputs

- aa_num -amino acid numbers, vector
- aa_typ -amino acid types, cell array of strings
- pdb_id -atom identifiers using PDB nomenclature,
- cell array of strings
- chems -chemical shifts in ppm
- Note: direct calls are discouraged, see the protein HOWTO docu-
- ment for instructions on importing protein data.

## Implementation structure

- Reads BMRB file, extracts amino acid numbers, amino acid types,
- PDB atom identifiers and chemical shifts. Syntax:
- [aa_num,aa_typ,pdb_id,chemsh]=read_bmrb(bmrb_file_name)
- bmrb_file_name -a character string with the file name
- aa_num -amino acid numbers, vector
- aa_typ -amino acid types, cell array of strings
- pdb_id -atom identifiers using PDB nomenclature,
- cell array of strings
- chems -chemical shifts in ppm
- Note: direct calls are discouraged, see the protein HOWTO docu-
- ment for instructions on importing protein data.
- Check consistency

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fopen()`, `feof()`, `fgetl()`, `regexprep()`, `textscan()`, `all()`, `cellfun()`, `aa_num()`, `chemsh()`, `upper()`, `fclose()`, `ischar()`.
