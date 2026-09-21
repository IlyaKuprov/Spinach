# interfaces/pdb_bmrb/protein.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/pdb_bmrb/protein.m`
- Signature: `[sys,inter,aux]=protein(pdb_file,bmrb_file,options)`
- Total lines: 524

## Purpose

Protein data import function. Parses PDB and BMRB data, runs a J-coupl- ing guess, a CSA guess and outputs Spinach data structures. Syntax: [sys,inter]=protein(pdb_file,bmrb_file,options)

## Physical / mathematical content

- PDB/BMRB interfaces. These files bridge biomolecular structure/assignment data and Spinach input structures, including atom selection, coordinates, and chemical-shift metadata.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Parameters / inputs

- pdb_file -string containing the name of the PDB file
- bmrb_file -string containing the name of the BMRB file
- options.select -'backbone' imports protein backbone up to
- CB and HB, 'backbone-minimal' only imports
- the backbone, 'backbone-hsqc' is the same
- as backbone, but with GLN and ASN side chain
- amide groups included, 'all' imports every-
- thing that is assigned in BMRB. If a list of
- numbers is supplied, atoms with those serial
- numbers in the PDB file are imported; every
- number must be present in the file, atoms of
- unsupported types (oxygen, sulphur, OH pro-
- tons) are dropped with a warning, and those
- without a BMRB assignment are kept or dele-
- ted according to options.noshift.
- options.pdb_mol -the number of molecule if there are multiple
- molecules in the pdb file
- options.noshift -'keep' places unassigned atoms between -1 and
- 0 ppm, 'delete' removes them from the system
- options.deuterate -a cell array of character strings, replaces
- protons with the specified PDB identifiers
- with deuterons; 'non-Me' deuterates every-
- thing except methyl groups
- options.nh_csa -peptide bond CSAs differ across literature,
- the following options are available:
- 'bax' for H:[6.00 0.00 -6.00], N:[-108.0 62.0 46.0] ppm
- 'tcb' for H:[7.00 0.00 -7.00], N:[-125.0 45.0 80.0] ppm
- 'pol' for H:[6.66 0.66 -7.33], N:[ -92.4 34.7 57.7] ppm
- the default is 'tcb'.

## Outputs

- sys.isotopes -Nspins x 1 cell array of strings
- sys.labels -Nspins x 1 cell array of strings containing
- standard IUPAC protein atom labels
- inter.coordinates -Nspins x 3 matrix, Angstrom.
- inter.zeeman.iso -Nspins x 1 cell array of numbers, ppm.
- Isotropic chemical shifts go here.
- inter.zeeman.matrix -Nspins x 1 cell array of 3x3 matrices, ppm.
- Chemical shift anisotropies go here.
- inter.coupling.scalar -Nspins x Nspins cell array of scalar coup-
- lings, all in Hz.
- aux.pdb_aa_num -pdb amino acid number for each spin
- aux.pdb_aa_typ -pdb amino acid type for each spin

## Implementation structure

- Protein data import function. Parses PDB and BMRB data, runs a J-coupl-
- ing guess, a CSA guess and outputs Spinach data structures. Syntax:
- [sys,inter]=protein(pdb_file,bmrb_file,options)
- pdb_file -string containing the name of the PDB file
- bmrb_file -string containing the name of the BMRB file
- options.select -'backbone' imports protein backbone up to
- CB and HB, 'backbone-minimal' only imports
- the backbone, 'backbone-hsqc' is the same
- as backbone, but with GLN and ASN side chain
- amide groups included, 'all' imports every-
- thing that is assigned in BMRB. If a list of
- numbers is supplied, spins with those num-

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `exist()`, `isfield()`, `grumble()`, `read_pdb_pro()`, `read_bmrb()`, `ismember()`, `pdb_aa_num()`, `pdb_atom_id()`, `pdb_aa_typ()`, `pdb_coords()`, `bmrb_atom_id()`, `bmrb_chemsh()`, `all()`, `strcmp()`, `bmrb_aa_typ()`, `num2str()`.
