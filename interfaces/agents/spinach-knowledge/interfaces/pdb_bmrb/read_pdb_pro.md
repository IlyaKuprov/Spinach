# interfaces/pdb_bmrb/read_pdb_pro.m

- Signature: `[aa_num,aa_typ,pdb_id,coords,pdb_ser]=read_pdb_pro(pdb_file_name,mod_id)`

## Purpose

Reads a PDB file and returns residue numbers, residue types, PDB atom identifiers, Cartesian coordinates, and PDB atom serial numbers.

## Parameters / inputs

- `pdb_file_name` — character string naming the PDB file.
- `mod_id` — model number to read; must be a finite positive integer. Use `1` for a file without `MODEL` records.

## Outputs

All outputs have matching entries for lines accepted by the parser:

- `aa_num` — residue number for each atom.
- `aa_typ` — cell array of residue identifiers (for example, `'TYR'`), converted to uppercase.
- `pdb_id` — cell array of PDB atom identifiers (for example, `'HE2'`).
- `coords` — cell array of three-element Cartesian coordinate vectors in angstroms.
- `pdb_ser` — PDB atom serial number for each atom.

## Implementation structure

The function scans for a matching `MODEL` record. If the file has no `MODEL` records and `mod_id` is `1`, it rewinds and reads the file as a single model. If a requested model is absent, the scan reaches end of file and no records are parsed; the function does not raise a separate not-found error. The parser then applies a space-delimited `ATOM`-record field layout, expecting the atom serial number, atom name, residue name, residue number, coordinates, and three further numeric fields; parsing stops at `ENDMDL` or end of file. Outputs are converted to column vectors, and the file is closed before return.