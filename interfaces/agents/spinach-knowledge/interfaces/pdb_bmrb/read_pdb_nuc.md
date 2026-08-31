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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 36-37: Check consistency; implemented by `grumble(pdb_file_name)`.
- Lines 39-40: Open the PDB file; implemented by `file_id=fopen(pdb_file_name,'r')`.
- Lines 42-43: Get the outputs started; implemented by `res_num=[]; res_typ={}`.
- Lines 46-47: Parse the PDB file; implemented by `while ~feof(file_id)`.
- Lines 60-61: Capitalise amino acid type specifications; implemented by `for n=1:numel(res_typ)`.
- Lines 65-66: Make outputs column vectors; implemented by `res_num=res_num'; res_typ=res_typ'`.
- Lines 69-70: Close the PDB file; implemented by `fclose(file_id)`.

### Control flow inferred from the code

- Line 47: `while` loop over `~feof(file_id)`.
- Line 49: conditional branch on `~isempty(data_line)`.
- Line 51: conditional branch on `all(~cellfun(@isempty,parsed_string))`.
- Line 61: `for` loop over `n=1:numel(res_typ)`.

### Key state/data transformations

- Lines 40: computes `file_id` using `file_id=fopen(pdb_file_name,'r')`.
- Lines 43: computes `res_num` using `res_num=[]; res_typ={}`.
- Lines 44: computes `pdb_id` using `pdb_id={}; coords={}`.
- Lines 48: computes `data_line` using `data_line=fgetl(file_id)`.
- Lines 50: computes `parsed_string` using `parsed_string=textscan(data_line,'ATOM %f %s %s %f %f %f %f %f %f %s','delimiter',' ','MultipleDelimsAsOne',1)`.
- Lines 52: computes `res_num(end+1)` using `res_num(end+1)=parsed_string{4}`.
- Lines 53: computes `res_typ{end+1}` using `res_typ{end+1}=parsed_string{3}{1}`.
- Lines 54: computes `pdb_id{end+1}` using `pdb_id{end+1}=parsed_string{2}{1}`.
- Lines 55: computes `coords{end+1}` using `coords{end+1}=[parsed_string{5:7}]`.
- Lines 62: computes `res_typ{n}` using `res_typ{n}=upper(res_typ{n})`.

### Local helper functions

- Line 75: `grumble()` — `function grumble(pdb_file_name)`. I think that it might well be premature to make an award of a Prize to Watson and Crick, because of existing uncertainty
  - Representative operation: `if ~ischar(pdb_file_name)`.
  - Representative operation: `error('pdb_file_name must be a character string.')`.

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
- one model.

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
