# interfaces/pdb_bmrb/nuclacid.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/pdb_bmrb/nuclacid.m`
- Signature: `[sys,inter]=nuclacid(pdb_file,shift_file,options)`
- Total lines: 224

## Purpose

Nucleic acid data import function. Parses PDB and chemical shift data, runs a J-coupling guess using guess_j_nuc.m function and outputs sys and inter data structures that are required by the create.m gateway function in Spinach. Syntax: [sys,inter]=nuclacid(pdb_file,shift_file,options)

## Physical / mathematical content

- PDB/BMRB interfaces. These files bridge biomolecular structure/assignment data and Spinach input structures, including atom selection, coordinates, and chemical-shift metadata.
- Chemical-shift anisotropy is present: shielding is treated as a second-rank tensor whose orientation relative to the field or rotor axis modulates line shapes and transfer dynamics.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `grumble()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 51-52: Check consistency; implemented by `grumble(pdb_file,shift_file,options)`.
- Lines 54-55: Parse the PDB file; implemented by `[pdb_res_num,pdb_res_typ,pdb_atom_id,pdb_coords]=read_pdb_nuc(pdb_file)`.
- Lines 57-58: Open the shift file; implemented by `file_id=fopen(shift_file,'r')`.
- Lines 60-61: Get the outputs started; implemented by `bmrb_res_num=[]; bmrb_atom_id={}; bmrb_chemsh={}`.
- Lines 63-64: Parse the shift file; implemented by `while ~feof(file_id)`.
- Lines 76-77: Make outputs column vectors; implemented by `bmrb_res_num=bmrb_res_num'`.
- Lines 81-82: Close the shift file; implemented by `fclose(file_id)`.
- Lines 84-85: Remove oxygen, phosphorus and the exchangeable protons; implemented by `kill_mask=ismember(pdb_atom_id,{'O5''','O4''','O4','O6','O2','O2''','O3''','O1P','O2P','P','H2''','H5T','HO''2'})`.
- Lines 90-91: Match chemical shifts; implemented by `pdb_chemsh=cell(numel(pdb_atom_id),1)`.
- Lines 95-96: Pull the current amino acid from BMRB; implemented by `select_mask=(bmrb_res_num==pdb_res_num(n))`.
- Lines 100-101: Match PDB and BMRB data; implemented by `if ismember(pdb_atom_id{n},bmrb_atoms)`.
- Lines 108-109: Replace primes with 'p' symbols for convenience; implemented by `for n=1:numel(pdb_atom_id)`.
- Lines 113-114: Estimate J-couplings; implemented by `scalar_couplings=guess_j_nuc(pdb_res_num,pdb_res_typ,pdb_atom_id,pdb_coords)`.
- Lines 116-119: CSA estimation goes here; implemented by `isotopes=cell(1,numel(pdb_atom_id))`.
- Lines 118-119: Assign isotopes and labels; implemented by `isotopes=cell(1,numel(pdb_atom_id))`.
- Lines 135-136: Find missing chemical shifts; implemented by `missing_shifts=find(cellfun(@isempty,pdb_chemsh))'`.
- Lines 139-140: Process unassigned chemical shifts; implemented by `subset=true(size(pdb_atom_id))`.
- Lines 143-144: Put unassigned spins between -1.0 and 0.0 ppm; implemented by `erzatz_shifts=linspace(-1,0,numel(missing_shifts))`.

### Control flow inferred from the code

- Line 64: `while` loop over `~feof(file_id)`.
- Line 66: conditional branch on `~isempty(data_line)`.
- Line 68: conditional branch on `all(~cellfun(@isempty,parsed_string))`.
- Line 93: `for` loop over `n=1:numel(pdb_atom_id)`.
- Line 101: conditional branch on `ismember(pdb_atom_id{n},bmrb_atoms)`.
- Line 109: `for` loop over `n=1:numel(pdb_atom_id)`.
- Line 120: `for` loop over `n=1:numel(pdb_atom_id)`.
- Line 121: dispatches on `pdb_atom_id{n}(1)`; cases `'H'`, `'C'`, `'N'`, `'P'`.
- Line 141: conditional branch on `strcmp(options.noshift,'keep')`.
- Line 145: `for` loop over `n=1:numel(missing_shifts)`.
- Line 154: `for` loop over `n=missing_shifts`.
- Line 178: `for` loop over `n=1:numel(pdb_atom_id)`.
- Line 184: `for` loop over `n=1:numel(sys.labels)`.
- Line 185: conditional branch on `ismember([pdb_res_typ{n} ':' pdb_atom_id{n}],options.deut_list)`.

### Key state/data transformations

- Lines 55: computes `[pdb_res_num,pdb_res_typ,pdb_atom_id,pdb_coords]` using `[pdb_res_num,pdb_res_typ,pdb_atom_id,pdb_coords]=read_pdb_nuc(pdb_file)`.
- Lines 58: computes `file_id` using `file_id=fopen(shift_file,'r')`.
- Lines 61: computes `bmrb_res_num` using `bmrb_res_num=[]; bmrb_atom_id={}; bmrb_chemsh={}`.
- Lines 65: computes `data_line` using `data_line=fgetl(file_id)`.
- Lines 67: computes `parsed_string` using `parsed_string=textscan(data_line,'%f %s %f','delimiter','\t','MultipleDelimsAsOne',1)`.
- Lines 69: computes `bmrb_res_num(end+1)` using `bmrb_res_num(end+1)=parsed_string{1}`.
- Lines 70: computes `bmrb_atom_id{end+1}` using `bmrb_atom_id{end+1}=parsed_string{2}{1}`.
- Lines 71: computes `bmrb_chemsh{end+1}` using `bmrb_chemsh{end+1}=parsed_string{3}`.
- Lines 78: computes `bmrb_atom_id` using `bmrb_atom_id=bmrb_atom_id'`.
- Lines 79: computes `bmrb_chemsh` using `bmrb_chemsh=bmrb_chemsh'`.
- Lines 85: computes `kill_mask` using `kill_mask=ismember(pdb_atom_id,{'O5''','O4''','O4','O6','O2','O2''','O3''','O1P','O2P','P','H2''','H5T','HO''2'})`.
- Lines 86: computes `pdb_res_num(kill_mask)` using `pdb_res_num(kill_mask)=[]; pdb_atom_id(kill_mask)=[]`.
- Lines 87: computes `pdb_res_typ(kill_mask)` using `pdb_res_typ(kill_mask)=[]; pdb_coords(kill_mask)=[]`.
- Lines 91: computes `pdb_chemsh` using `pdb_chemsh=cell(numel(pdb_atom_id),1)`.
- Lines 97: computes `bmrb_atoms` using `bmrb_atoms=bmrb_atom_id(select_mask)`.
- Lines 98: computes `bmrb_shifts` using `bmrb_shifts=bmrb_chemsh(select_mask)`.
- Lines 102: computes `pdb_chemsh{n}` using `pdb_chemsh{n}=bmrb_shifts(strcmp(pdb_atom_id{n},bmrb_atoms)); pdb_chemsh{n}=pdb_chemsh{n}{1}`.
- Lines 110: computes `pdb_atom_id{n}` using `pdb_atom_id{n}=strrep(pdb_atom_id{n},'''','p')`.

### Local helper functions

- Line 200: `grumble()` — `function grumble(pdb_file,shift_file,options)`.
  - Representative operation: `if ~ischar(pdb_file)`.
  - Representative operation: `error('pdb_file must be a character string specifying a file name.')`.

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

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `read_pdb_nuc()`, `fopen()`, `feof()`, `fgetl()`, `textscan()`, `all()`, `cellfun()`, `bmrb_res_num()`, `fclose()`, `ismember()`, `pdb_res_num()`, `pdb_atom_id()`, `pdb_res_typ()`, `pdb_coords()`, `bmrb_atom_id()`.
