# interfaces/pdb_bmrb/read_pdb_pro.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/pdb_bmrb/read_pdb_pro.m`
- Signature: `[aa_num,aa_typ,pdb_id,coords]=read_pdb_pro(pdb_file_name,mod_id)`
- Total lines: 109

## Purpose

Reads the a PDB file and returns amino acid numbers, amino acid types, PDB atom identifiers and Cartesian coordinates. Syntax: [aa_num,aa_typ,pdb_id,coords]=read_pdb_pro(pdb_file_name,instance)

## Physical / mathematical content

- PDB/BMRB interfaces. These files bridge biomolecular structure/assignment data and Spinach input structures, including atom selection, coordinates, and chemical-shift metadata.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 35-36: Check consistency; implemented by `grumble(pdb_file_name,mod_id)`.
- Lines 38-39: Open the PDB file; implemented by `file_id=fopen(pdb_file_name,'r')`.
- Lines 41-42: Scroll to the selected structure; implemented by `while ~feof(file_id)`.
- Lines 49-50: Get the outputs started; implemented by `aa_num=[]; aa_typ={}`.
- Lines 53-54: Parse the PDB file; implemented by `while ~feof(file_id)`.
- Lines 68-69: Capitalize amino acid type specifications; implemented by `for n=1:numel(aa_typ)`.
- Lines 73-74: Make outputs column vectors; implemented by `aa_num=aa_num'; aa_typ=aa_typ'`.
- Lines 77-78: Close the PDB file; implemented by `fclose(file_id)`.

### Control flow inferred from the code

- Line 42: `while` loop over `~feof(file_id)`.
- Line 44: conditional branch on `parsed_string{1}==mod_id`.
- Line 54: `while` loop over `~feof(file_id)`.
- Line 56: conditional branch on `(numel(data_line)>=6)&&strcmp('ENDMDL',data_line(1:6)), break; end`.
- Line 57: conditional branch on `~isempty(data_line)`.
- Line 59: conditional branch on `all(~cellfun(@isempty,parsed_string))`.
- Line 69: `for` loop over `n=1:numel(aa_typ)`.

### Key state/data transformations

- Lines 39: computes `file_id` using `file_id=fopen(pdb_file_name,'r')`.
- Lines 43: computes `parsed_string` using `parsed_string=textscan(fgetl(file_id),'MODEL %d','delimiter',' ','MultipleDelimsAsOne',1)`.
- Lines 50: computes `aa_num` using `aa_num=[]; aa_typ={}`.
- Lines 51: computes `pdb_id` using `pdb_id={}; coords={}`.
- Lines 55: computes `data_line` using `data_line=fgetl(file_id)`.
- Lines 60: computes `aa_num(end+1)` using `aa_num(end+1)=parsed_string{5}`.
- Lines 61: computes `aa_typ{end+1}` using `aa_typ{end+1}=parsed_string{3}{1}`.
- Lines 62: computes `pdb_id{end+1}` using `pdb_id{end+1}=parsed_string{2}{1}`.
- Lines 63: computes `coords{end+1}` using `coords{end+1}=[parsed_string{6:8}]`.
- Lines 70: computes `aa_typ{n}` using `aa_typ{n}=upper(aa_typ{n})`.

### Local helper functions

- Line 83: `grumble()` — `function grumble(pdb_file_name,mod_id)`. Heartless sterility, obliteration of all melody, all tonal charm, all
  - Representative operation: `if ~ischar(pdb_file_name)`.
  - Representative operation: `error('pdb_file_name must be a character string.')`.

## Parameters / inputs

- pdb_file_name -a character string with the file name
- mod_id -the number of model that should be
- read from the pdb file

## Outputs

- aa_num -nspins x 1 vector giving the number of
- the amino acid to which each spin belongs
- aa_typ -nspins x 1 cell array of strings giving
- the PDB identifier of the amino acid to
- which each spin belongs (e.g. 'TYR')
- pdb_id -nspins x 1 cell array of strings giving
- the PDB identifier of the protein atom
- type to which each spin belongs (e.g. 'HE2')
- coords -nspins x 1 cell array of 3-vectors giving
- Cartesian coordinates of each spin in Angstrom

## Implementation structure

- Reads the a PDB file and returns amino acid numbers, amino acid types,
- PDB atom identifiers and Cartesian coordinates. Syntax:
- [aa_num,aa_typ,pdb_id,coords]=read_pdb_pro(pdb_file_name,instance)
- pdb_file_name -a character string with the file name
- mod_id -the number of model that should be
- read from the pdb file
- aa_num -nspins x 1 vector giving the number of
- the amino acid to which each spin belongs
- aa_typ -nspins x 1 cell array of strings giving
- the PDB identifier of the amino acid to
- which each spin belongs (e.g. 'TYR')
- pdb_id -nspins x 1 cell array of strings giving

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `fopen()`, `feof()`, `textscan()`, `fgetl()`, `num2str()`, `strcmp()`, `data_line()`, `all()`, `cellfun()`, `aa_num()`, `upper()`, `fclose()`, `ischar()`, `isscalar()`.
