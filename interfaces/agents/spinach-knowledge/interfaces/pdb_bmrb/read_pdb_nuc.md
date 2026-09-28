# interfaces/pdb_bmrb/read_pdb_nuc.m

- Signature: `[res_num,res_typ,pdb_id,coords]=read_pdb_nuc(pdb_file_name)`

## Purpose

Reads atom records from a PDB file and returns each parsed atom's residue number, residue type, PDB atom label, and Cartesian coordinates.

## Parameters / inputs

- `pdb_file_name` — a character string specifying the PDB file. Other input types are rejected.

## Outputs

- `res_num` — `nspins x 1` vector of nucleotide residue numbers.
- `res_typ` — `nspins x 1` cell array of nucleotide PDB identifiers (for example, `'GUA'`), converted to uppercase.
- `pdb_id` — `nspins x 1` cell array of nucleic-acid atom PDB identifiers (for example, `'C1P'`).
- `coords` — `nspins x 1` cell array of three-element Cartesian coordinate vectors in angstroms.

## Parsing and validation

The parser reads lines matching its `ATOM` field pattern, taking the residue number, residue type, atom label, and three coordinate values from each successfully parsed line. It removes the chain-identifier column before parsing, accepts a chain identifier without returning it, and rejects files whose parsed atoms have more than one distinct chain identifier because residue numbers would be ambiguous. Prepare a PDB file containing one model and one chain; the code checks for multiple chains but does not separately validate the number of models. Outputs are returned as column vectors.