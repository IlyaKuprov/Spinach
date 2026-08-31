# etc/estimators/guess_j_nuc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/estimators/guess_j_nuc.m`
- Signature: `jmatrix=guess_j_nuc(nuc_num,nuc_typ,pdb_id,coords)`
- Total lines: 462

## Purpose

RNA assignments of J-couplings from literature values and Karplus cur- ves. Syntax: jmatrix=guess_j_nuc(nuc_num,nuc_typ,pdb_id,coords)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 37-38: Check consistency; implemented by `grumble(nuc_num,nuc_typ,pdb_id,coords)`.
- Lines 40-41: Preallocate the answer; implemented by `jmatrix=cell(numel(coords),numel(coords))`.
- Lines 43-44: Number the atoms; implemented by `numbers=1:numel(coords)`.
- Lines 46-47: Determine the connectivity graph; implemented by `proxmatrix=false(numel(coords),numel(coords))`.
- Lines 54-55: Get all connected subgraphs up to size 2; implemented by `subgraphs=dfpt(sparse(proxmatrix),2)`.
- Lines 57-58: Set generic coupling values; implemented by `J_NH=86.0; J_CN=20.0; J_CH=180.0; J_CC=65.0`.
- Lines 60-61: Spec all ordered connected pairs; implemented by `pairs_database={`.
- Lines 63-64: RNA specific values from the literature; implemented by `'C1p_N9' , 11.0; 'C1p_N1' , 12.0; 'H_N1' , 95.0; 'H3_N3' , 91.0`.
- Lines 70-71: Generics for the conjugated ring (to be replaced with DFT values -Zenawi); implemented by `'C2_H2' , J_CH; 'C5_N7' , J_CN; 'C8_N7' , J_CN; 'C6_N6' , J_CN`.
- Lines 75-76: Generics for the sugar ring (to be replaced with DFT values -Zenawi); implemented by `'C1p_C2p' , J_CC; 'C2p_H2p' , J_CH; 'C3p_C4p' , J_CC; 'C4p_H4p' , J_CH`.
- Lines 81-82: Loop over subgraphs; implemented by `disp(' '); disp('############# ONE-BOND J-COUPLING SUMMARY ###############')`.
- Lines 86-87: Extract labels, residues and numbers; implemented by `spin_labels=pdb_id(subgraphs(n,:))`.
- Lines 92-93: Index the spins; implemented by `[spin_labels,index]=sort(spin_labels)`.
- Lines 98-99: Consult the database; implemented by `if isscalar(spin_labels)`.
- Lines 101-102: Inform the user about solitary spins; implemented by `disp(['WARNING: solitary spin detected, ' spin_resnames{1} '(' num2str(spin_resnums(1)) '):' spin_labels{1}])`.
- Lines 106-107: Build specification string; implemented by `spec_string=pairs_database(strcmp([spin_labels{1} '_' spin_labels{2}],pairs_database(:,1)),:)`.
- Lines 109-110: Assign the coupling; implemented by `if isempty(spec_string)`.
- Lines 112-114: Complain if nothing found; implemented by `disp(['WARNING: unknown atom pair or unphysical proximity, ' spin_labels{1} '_' spin_labels{2} ' in residues ' spin_resnames{1} '(' num2str(spin_resnums(1)) '), ' spin_r…`.

### Control flow inferred from the code

- Line 48: `for` loop over `n=1:numel(coords)`.
- Line 49: `for` loop over `k=1:numel(coords)`.
- Line 50: conditional branch on `(norm(coords{n}-coords{k},2)<1.55), proxmatrix(n,k)=1; end`.
- Line 84: `for` loop over `n=1:size(subgraphs,1)`.
- Line 99: conditional branch on `isscalar(spin_labels)`.
- Line 110: conditional branch on `isempty(spec_string)`.
- Line 193: `for` loop over `n=1:size(subgraphs,1)`.
- Line 208: conditional branch on `isscalar(spin_labels)`.
- Line 224: conditional branch on `isempty(spec_string)`.
- Line 243: conditional branch on `proxmatrix(spin_a,spin_b)&&proxmatrix(spin_b,spin_c)`.
- Line 246: conditional branch on `(~isempty(jmatrix{spin_a,spin_c}))&&(~ismember([spin_labels{1} '_' spin_labels{2} '_' spin_labels{3}],allowed_collision…`.
- Line 279: `for` loop over `n=1:size(subgraphs,1)`.
- Line 280: conditional branch on `any(sum(proxmatrix(subgraphs(n,:),subgraphs(n,:)))==4)`.
- Line 344: `for` loop over `n=1:size(subgraphs,1)`.

### Key state/data transformations

- Lines 41: computes `jmatrix` using `jmatrix=cell(numel(coords),numel(coords))`.
- Lines 44: computes `numbers` using `numbers=1:numel(coords)`.
- Lines 47: computes `proxmatrix` using `proxmatrix=false(numel(coords),numel(coords))`.
- Lines 55: computes `subgraphs` using `subgraphs=dfpt(sparse(proxmatrix),2)`.
- Lines 58: computes `J_NH` using `J_NH=86.0; J_CN=20.0; J_CH=180.0; J_CC=65.0`.
- Lines 61: computes `pairs_database` using `pairs_database={`.
- Lines 87: computes `spin_labels` using `spin_labels=pdb_id(subgraphs(n,:))`.
- Lines 88: computes `spin_numbers` using `spin_numbers=numbers(subgraphs(n,:))`.
- Lines 89: computes `spin_resnames` using `spin_resnames=nuc_typ(subgraphs(n,:))`.
- Lines 90: computes `spin_resnums` using `spin_resnums=nuc_num(subgraphs(n,:))`.
- Lines 93: computes `[spin_labels,index]` using `[spin_labels,index]=sort(spin_labels)`.
- Lines 107: computes `spec_string` using `spec_string=pairs_database(strcmp([spin_labels{1} '_' spin_labels{2}],pairs_database(:,1)),:)`.
- Lines 124: computes `jmatrix{spin_numbers(1),spin_numbers(2)}` using `jmatrix{spin_numbers(1),spin_numbers(2)}=spec_string{2}`.
- Lines 141: computes `J_CCC` using `J_CCC=-1.2`.
- Lines 142: computes `J_CCH` using `J_CCH=0`.
- Lines 143: computes `J_CCN` using `J_CCN=7.0`.
- Lines 144: computes `J_NCN` using `J_NCN=0`.
- Lines 145: computes `J_NCH` using `J_NCH=0`.

### Local helper functions

- Line 436: `grumble()` — `function grumble(nuc_num,nuc_typ,pdb_id,coords)`.
  - Representative operation: `if ~isnumeric(nuc_num)`.
  - Representative operation: `error('nuc_num must be a vector of positive integers.')`.

## Parameters / inputs

- nuc_num -a vector of nucleotide numbers
- nuc_typ -a cell array of nucleotide types
- pdb_id -a cell array of PDB atom identifiers
- coords -a cell array of coordinate vectors

## Outputs

- jmatrix -a cell array of J-couplings in Hz
- Database nomenclature is as follows:
- 1. Atoms in the subgraph descriptor are listed alphabetically to make
- the descriptors unique.
- 2. The four numbers refer to the bonding order, e.g. 1, 3, 2, 4 means
- that the first atom in the descriptor is bonded to the third, which
- is bonded to the second, which is bonded to the fourth. The coupling
- in this case is between atom 1 and atom 4 in the descriptor.

## Implementation structure

- RNA assignments of J-couplings from literature values and Karplus cur-
- ves. Syntax:
- jmatrix=guess_j_nuc(nuc_num,nuc_typ,pdb_id,coords)
- nuc_num -a vector of nucleotide numbers
- nuc_typ -a cell array of nucleotide types
- pdb_id -a cell array of PDB atom identifiers
- coords -a cell array of coordinate vectors
- jmatrix -a cell array of J-couplings in Hz
- Database nomenclature is as follows:
- 1. Atoms in the subgraph descriptor are listed alphabetically to make
- the descriptors unique.
- 2. The four numbers refer to the bonding order, e.g. 1, 3, 2, 4 means

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `false()`, `proxmatrix()`, `dfpt()`, `pdb_id()`, `subgraphs()`, `numbers()`, `nuc_typ()`, `nuc_num()`, `spin_numbers()`, `spin_resnames()`, `spin_resnums()`, `isscalar()`, `num2str()`, `pairs_database()`, `strcmp()`.
