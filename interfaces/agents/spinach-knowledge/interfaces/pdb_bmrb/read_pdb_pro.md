# interfaces/pdb_bmrb/read_pdb_pro.m

- Signature: `[aa_num,aa_typ,pdb_id,coords,pdb_ser]=read_pdb_pro(pdb_file_name,mod_id)`

## Purpose

Reads the a PDB file and returns amino acid numbers, amino acid types, PDB atom identifiers and Cartesian coordinates. Syntax: [aa_num,aa_typ,pdb_id,coords,pdb_ser]=read_pdb_pro(pdb_file_name,mod_id)

## Physical / mathematical content

- PDB/BMRB interfaces. These files bridge biomolecular structure/assignment data and Spinach input structures, including atom selection, coordinates, and chemical-shift metadata.

## Numerical / algorithmic content

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
- pdb_ser -nspins x 1 vector giving the PDB atom serial
- number of each spin

## Implementation structure

- Reads the a PDB file and returns amino acid numbers, amino acid types,
- PDB atom identifiers and Cartesian coordinates. Syntax:
- [aa_num,aa_typ,pdb_id,coords,pdb_ser]=read_pdb_pro(pdb_file_name,mod_id)
- pdb_file_name -a character string with the file name
- mod_id -the number of model that should be
- read from the pdb file
- aa_num -nspins x 1 vector giving the number of
- the amino acid to which each spin belongs
- aa_typ -nspins x 1 cell array of strings giving
- the PDB identifier of the amino acid to
- which each spin belongs (e.g. 'TYR')
- pdb_id -nspins x 1 cell array of strings giving
