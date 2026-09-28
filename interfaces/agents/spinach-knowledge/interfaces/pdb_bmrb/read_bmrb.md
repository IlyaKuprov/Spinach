# interfaces/pdb_bmrb/read_bmrb.m

- Signature: `[aa_num,aa_typ,pdb_id,chemsh]=read_bmrb(bmrb_file_name)`

## Purpose

Reads a BMRB file and extracts amino acid numbers, amino acid types, atom identifiers using PDB nomenclature, and chemical shifts.

## Physical / mathematical content

- PDB/BMRB interfaces bridge biomolecular structure and assignment data with Spinach input structures, including atom selection, coordinates, and chemical-shift metadata.

## Numerical / algorithmic content

- Reads the file line by line, replaces tabs with spaces, and skips empty lines.
- Parses space-delimited lines as eight fields with types `%f %f %s %s %s %f %f %f`, treating consecutive delimiters as one. A line contributes an output record only if all eight parsed fields are nonempty.
- Uses fields 2, 3, 4, and 6 for `aa_num`, `aa_typ`, `pdb_id`, and `chemsh`, respectively. Converts amino acid types to uppercase and returns all four outputs as column vectors.
- Raises an error if no chemical shifts were parsed, reporting that the data format is not canonical.

## Parameters / inputs

- `bmrb_file_name` — file name as a character string. Other input types are rejected.

## Outputs

- `aa_num` — amino acid numbers, numeric column vector.
- `aa_typ` — amino acid types, column cell array of uppercase strings.
- `pdb_id` — atom identifiers using PDB nomenclature, column cell array of strings.
- `chemsh` — chemical shifts in ppm, numeric column vector.

## Implementation structure

- Validates the input type, opens the BMRB file for reading, initializes the outputs, parses eligible lines, converts amino acid types to uppercase, transposes the outputs into column vectors, and closes the file before checking for an empty result.
- Direct calls are discouraged; see the protein HOWTO document for instructions on importing protein data.
