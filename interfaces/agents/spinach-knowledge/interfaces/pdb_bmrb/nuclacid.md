# interfaces/pdb_bmrb/nuclacid.m

- Signature: `[sys,inter]=nuclacid(pdb_file,shift_file,options)`

## Purpose

Nucleic acid data import function. Parses PDB and chemical shift data, runs a J-coupling guess using guess_j_nuc.m function and outputs sys and inter data structures that are required by the create.m gateway function in Spinach. Syntax: [sys,inter]=nuclacid(pdb_file,shift_file,options)

## Physical / mathematical content

- PDB/BMRB interfaces. These files bridge biomolecular structure/assignment data and Spinach input structures, including atom selection, coordinates, and chemical-shift metadata.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

## Parameters / inputs

- pdb_file -a character string containing the name
- of the PDB file
- shift_file -a character string containing the name
- of the chemical shift file, ASCII for-
- matted as [residue_number atom_id shift],
- see example.txt in examples/nmr_nucleic
- options.deut_list -a cell array of strings, specifying which
- atoms should be assumed to be deuterated,
- for example {'ADE:H2pp'}. When an atom is
- deuterated, J-couplings are reduced appro-
- priately.
- options.noshift -'keep' places unassigned atoms between -1
- and 0 ppm, 'delete' removes them from the
- system

## Outputs

- sys.isotopes -Nspins x 1 cell array of strings
- sys.labels -Nspins x 1 cell array of strings
- containing standard IUPAC DNA/RNA
- atom labels
- inter.coordinates -Nspins x 3 matrix, Angstrom.
- inter.zeeman.scalar -Nspins x 1 cell array of numbers,
- ppm. Isotropic chemical shifts go
- here.
- inter.coupling.scalar -Nspins x Nspins cell array of sca-
- lar couplings, all in Hz.

## Implementation structure

- Nucleic acid data import function. Parses PDB and chemical shift
- data, runs a J-coupling guess using guess_j_nuc.m function and
- outputs sys and inter data structures that are required by the
- create.m gateway function in Spinach. Syntax:
- [sys,inter]=nuclacid(pdb_file,shift_file,options)
- pdb_file -a character string containing the name
- of the PDB file
- shift_file -a character string containing the name
- of the chemical shift file, ASCII for-
- matted as [residue_number atom_id shift],
- see example.txt in examples/nmr_nucleic
- options.deut_list -a cell array of strings, specifying which
