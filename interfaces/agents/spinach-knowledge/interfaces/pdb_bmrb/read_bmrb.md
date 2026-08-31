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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 30-31: Check consistency; implemented by `grumble(bmrb_file_name)`.
- Lines 33-34: Open the BMRB file; implemented by `file_id=fopen(bmrb_file_name,'r')`.
- Lines 36-37: Get the outputs started; implemented by `aa_num=[]; aa_typ={}`.
- Lines 40-41: Parse the BMRB file; implemented by `while ~feof(file_id)`.
- Lines 43-44: Get a line; implemented by `data_line=fgetl(file_id)`.
- Lines 46-47: Replace tabs with spaces; implemented by `data_line=regexprep(data_line,'\t',' ')`.
- Lines 49-50: If not empty, read; implemented by `if ~isempty(data_line)`.
- Lines 62-63: Capitalize amino acid type specifications; implemented by `for n=1:numel(aa_typ)`.
- Lines 67-68: Make outputs column vectors; implemented by `aa_num=aa_num'; aa_typ=aa_typ'`.
- Lines 71-72: Close the BMRB file; implemented by `fclose(file_id)`.
- Lines 74-75: Check for empty returns; implemented by `if isempty(chemsh)`.

### Control flow inferred from the code

- Line 41: `while` loop over `~feof(file_id)`.
- Line 50: conditional branch on `~isempty(data_line)`.
- Line 52: conditional branch on `all(~cellfun(@isempty,parsed_string))`.
- Line 63: `for` loop over `n=1:numel(aa_typ)`.
- Line 75: conditional branch on `isempty(chemsh)`.

### Key state/data transformations

- Lines 34: computes `file_id` using `file_id=fopen(bmrb_file_name,'r')`.
- Lines 37: computes `aa_num` using `aa_num=[]; aa_typ={}`.
- Lines 38: computes `pdb_id` using `pdb_id={}; chemsh=[]`.
- Lines 44: computes `data_line` using `data_line=fgetl(file_id)`.
- Lines 51: computes `parsed_string` using `parsed_string=textscan(data_line,'%f %f %s %s %s %f %f %f','delimiter',' ','MultipleDelimsAsOne',1)`.
- Lines 53: computes `aa_num(end+1)` using `aa_num(end+1)=parsed_string{2}`.
- Lines 54: computes `aa_typ{end+1}` using `aa_typ{end+1}=parsed_string{3}{1}`.
- Lines 55: computes `pdb_id{end+1}` using `pdb_id{end+1}=parsed_string{4}{1}`.
- Lines 56: computes `chemsh(end+1)` using `chemsh(end+1)=parsed_string{6}`.
- Lines 64: computes `aa_typ{n}` using `aa_typ{n}=upper(aa_typ{n})`.

### Local helper functions

- Line 82: `grumble()` — `function grumble(bmrb_file_name)`. If it's not true, it's well invented. Dante
  - Representative operation: `if ~ischar(bmrb_file_name)`.
  - Representative operation: `error('bmrb_file_name must be a character string.')`.

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
