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

## Code-derived implementation details

### Comment-guided execution stages

- Lines 78-79: Set defaults; implemented by `if ~exist('options','var')`.
- Lines 89-90: Check consistency; implemented by `grumble(pdb_file,bmrb_file,options)`.
- Lines 92-93: Parse the PDB file; implemented by `[pdb_aa_num,pdb_aa_typ,pdb_atom_id,pdb_coords,pdb_ser]=read_pdb_pro(pdb_file,options.pdb_mol)`.
- Lines 103-104: Parse the BMRB file; implemented by `[bmrb_aa_num,bmrb_aa_typ,bmrb_atom_id,bmrb_chemsh]=read_bmrb(bmrb_file)`.
- Lines 106-108: Remove oxygens, sulphurs and terminal atoms; implemented by `kill_mask=ismember(pdb_atom_id,{'O','OE','OE1','OE2','OD1','OD2','OG','OG1','HG1', 'OG2','OH','HH','SD','SG','OXT','O''','O'''''})`.
- Lines 121-122: Match chemical shifts; implemented by `pdb_chemsh=cell(numel(pdb_atom_id),1)`.
- Lines 126-127: Pull the current amino acid from BMRB; implemented by `select_mask=(bmrb_aa_num==pdb_aa_num(n))`.
- Lines 131-132: Make sure amino acid types match; implemented by `if ~all(strcmp(pdb_aa_typ{n},bmrb_aa_typ(select_mask)))`.
- Lines 139-140: Ugly heuristics (sorry!) to match PDB and BMRB data; implemented by `if ismember(pdb_atom_id{n},bmrb_atoms)`.
- Lines 142-143: If the atom is found, just import the shift; implemented by `pdb_chemsh{n}=bmrb_shifts(strcmp(pdb_atom_id{n},bmrb_atoms))`.
- Lines 149-150: Ignore the N-terminal amino group nitrogen (usually hard to assign); implemented by `disp(['WARNING: N-terminal NH2 group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.'])`.
- Lines 154-155: Ignore OH protons in serine (usually rapid exchange); implemented by `disp(['WARNING: OH group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.'])`.
- Lines 159-160: Ignore OH protons in threonine (usually rapid exchange); implemented by `disp(['WARNING: OH group of ' pdb_aa_typ{n} '(' num2str(pdb_aa_num(n)) ')' ' ignored.'])`.
- Lines 165-166: Use the shift reported for one methyl group proton for all three; implemented by `pdb_chemsh{n}=bmrb_shifts(strcmp('HB2',bmrb_atoms))`.
- Lines 173-174: Use the shift reported for one methyl group proton for all three; implemented by `pdb_chemsh{n}=bmrb_shifts(strcmp('HB',bmrb_atoms))`.
- Lines 181-182: Use the shift reported for one methyl group proton for all three; implemented by `pdb_chemsh{n}=bmrb_shifts(strcmp('HG2',bmrb_atoms))`.
- Lines 189-190: Use the shift reported for one methyl group proton for all three; implemented by `pdb_chemsh{n}=bmrb_shifts(strcmp('HG1',bmrb_atoms))`.
- Lines 197-198: Use the shift reported for one methyl group proton for all three; implemented by `pdb_chemsh{n}=bmrb_shifts(strcmp('HG12',bmrb_atoms))`.

### Control flow inferred from the code

- Line 75: conditional branch on `~exist('options','var')`.
- Line 81: conditional branch on `~isfield(options,'deuterate')`.
- Line 104: `for` loop over `n=1:numel(pdb_atom_id)`.
- Line 112: conditional branch on `~all(strcmp(pdb_aa_typ{n},bmrb_aa_typ(select_mask)))`.
- Line 120: conditional branch on `ismember(pdb_atom_id{n},bmrb_atoms)`.
- Line 287: `for` loop over `n=1:numel(pdb_atom_id)`.
- Line 288: dispatches on `pdb_atom_id{n}(1)`; cases `'H'`, `'C'`, `'N'`.
- Line 301: conditional branch on `isnumeric(options.select)`.
- Line 341: conditional branch on `strcmp(options.noshift,'keep')`.
- Line 345: `for` loop over `n=1:numel(missing_shifts)`.
- Line 356: `for` loop over `n=missing_shifts`.
- Line 371: conditional branch on `~isempty(options.deuterate)`.
- Line 377: conditional branch on `iscell(options.deuterate)`.
- Line 380: `for` loop over `n=1:numel(isotopes)`.

### Key state/data transformations

- Lines 80: computes `options.select` using `options.select='all'`.
- Lines 81: computes `options.pdb_mol` using `options.pdb_mol=1`.
- Lines 82: computes `options.noshift` using `options.noshift='keep'`.
- Lines 83: computes `options.deuterate` using `options.deuterate={}`.
- Lines 93: computes `[pdb_aa_num,pdb_aa_typ,pdb_atom_id,pdb_coords,pdb_ser]` using `[pdb_aa_num,pdb_aa_typ,pdb_atom_id,pdb_coords,pdb_ser]=read_pdb_pro(pdb_file,options.pdb_mol)`.
- Lines 104: computes `[bmrb_aa_num,bmrb_aa_typ,bmrb_atom_id,bmrb_chemsh]` using `[bmrb_aa_num,bmrb_aa_typ,bmrb_atom_id,bmrb_chemsh]=read_bmrb(bmrb_file)`.
- Lines 107-108: computes `kill_mask` using `kill_mask=ismember(pdb_atom_id,{'O','OE','OE1','OE2','OD1','OD2','OG','OG1','HG1', 'OG2','OH','HH','SD','SG','OXT','O''','O'''''})`.
- Lines 109: computes `pdb_aa_num(kill_mask)` using `pdb_aa_num(kill_mask)=[]; pdb_atom_id(kill_mask)=[]`.
- Lines 110: computes `pdb_aa_typ(kill_mask)` using `pdb_aa_typ(kill_mask)=[]; pdb_coords(kill_mask)=[]`.
- Lines 122: computes `pdb_chemsh` using `pdb_chemsh=cell(numel(pdb_atom_id),1)`.
- Lines 128: computes `bmrb_atoms` using `bmrb_atoms=bmrb_atom_id(select_mask)`.
- Lines 129: computes `bmrb_shifts` using `bmrb_shifts=bmrb_chemsh(select_mask)`.
- Lines 143: computes `pdb_chemsh{n}` using `pdb_chemsh{n}=bmrb_shifts(strcmp(pdb_atom_id{n},bmrb_atoms))`.
- Lines 287: computes `scalar_couplings` using `scalar_couplings=guess_j_pro(pdb_aa_num,pdb_aa_typ,pdb_atom_id,pdb_coords)`.
- Lines 290: computes `CSAs` using `CSAs=guess_csa_pro(pdb_aa_num,pdb_atom_id,pdb_coords,options)`.
- Lines 293: computes `isotopes` using `isotopes=cell(1,numel(pdb_atom_id))`.
- Lines 297: computes `isotopes{n}` using `isotopes{n}='1H'`.
- Lines 311: computes `subset` using `subset=ismember(pdb_ser,options.select)`, matching PDB atom serial numbers rather than positional indices; a requested serial number absent from the unfiltered PDB model raises an error earlier (lines 95-101), and a serial that is present yet names an atom type excluded by the oxygen/sulphur/terminal-atom kill mask (lines 106-112) is dropped from the subset with a printed warning (lines 113-119), not an error.

### Local helper functions

- Line 474: `grumble()` — `function grumble(pdb_file,bmrb_file,options)`.
  - Representative operation: `if ~ischar(pdb_file)`.
  - Representative operation: `error('pdb_file must be a character string specifying a file name.')`.

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
