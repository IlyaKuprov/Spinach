# etc/estimators/guess_j_pro.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/estimators/guess_j_pro.m`
- Signature: `jmatrix=guess_j_pro(aa_num,aa_typ,pdb_id,coords)`
- Total lines: 660

## Purpose

Assigns J-couplings from literature values and Karplus curves. Syntax: jmatrix=guess_j_pro(aa_num,aa_typ,pdb_id,coords)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

## Code-derived implementation details

### Comment-guided execution stages

- Lines 43-44: Check consistency; implemented by `grumble(aa_num,aa_typ,pdb_id,coords)`.
- Lines 46-47: Preallocate the answer; implemented by `jmatrix=cell(numel(coords),numel(coords))`.
- Lines 49-50: Number the atoms; implemented by `numbers=1:numel(coords)`.
- Lines 52-53: Determine the connectivity graph; implemented by `proxmatrix=false(numel(coords),numel(coords))`.
- Lines 60-61: Get all connected subgraphs up to size 2; implemented by `subgraphs=dfpt(sparse(proxmatrix),2)`.
- Lines 63-64: Set generic coupling values; implemented by `J_NH=-90.0; J_CN=-15.0; J_CH=140.0; J_CC=35.0`.
- Lines 66-67: Spec all ordered connected pairs; implemented by `pairs_database={`.
- Lines 69-70: Backbone has specific numbers; implemented by `'CA_CB' , +34.9; 'CA_HA' , +143.5; 'CA_N' , -10.7`.
- Lines 73-74: The rest is generic; implemented by `'CA_HA2' , J_CH; 'CA_HA3' , J_CH; 'CB_CG' , J_CC; 'CB_CG1' , J_CC`.
- Lines 93-94: Loop over subgraphs; implemented by `disp(' '); disp('############# ONE-BOND J-COUPLING SUMMARY ###############')`.
- Lines 97-98: Extract labels, residues and numbers; implemented by `spin_labels=pdb_id(subgraphs(n,:))`.
- Lines 103-104: Index the spins; implemented by `[spin_labels,index]=sort(spin_labels)`.
- Lines 109-110: Consult the database; implemented by `if isscalar(spin_labels)`.
- Lines 112-113: Inform the user about solitary spins; implemented by `disp(['WARNING: solitary spin detected, ' spin_resnames{1} '(' num2str(spin_resnums(1)) '):' spin_labels{1}])`.
- Lines 117-118: Build specification string; implemented by `spec_string=pairs_database(strcmp([spin_labels{1} '_' spin_labels{2}],pairs_database(:,1)),:)`.
- Lines 120-121: Assign the coupling; implemented by `if isempty(spec_string)`.
- Lines 123-125: Complain if nothing found; implemented by `disp(['WARNING: unknown atom pair or unphysical proximity, ' spin_labels{1} '_' spin_labels{2} ' in residues ' spin_resnames{1} '(' num2str(spin_resnums(1)) '), ' spin_r…`.
- Lines 129-130: Complain if multiple records found; implemented by `error([spec_string{1} ' coupling has been specified multiple times in the coupling database.'])`.

### Control flow inferred from the code

- Line 54: `for` loop over `n=1:numel(coords)`.
- Line 55: `for` loop over `k=1:numel(coords)`.
- Line 56: conditional branch on `(norm(coords{n}-coords{k},2)<1.60), proxmatrix(n,k)=1; end`.
- Line 95: `for` loop over `n=1:size(subgraphs,1)`.
- Line 110: conditional branch on `isscalar(spin_labels)`.
- Line 121: conditional branch on `isempty(spec_string)`.
- Line 220: `for` loop over `n=1:size(subgraphs,1)`.
- Line 235: conditional branch on `isscalar(spin_labels)`.
- Line 251: conditional branch on `isempty(spec_string)`.
- Line 270: conditional branch on `proxmatrix(spin_a,spin_b)&&proxmatrix(spin_b,spin_c)`.
- Line 273: conditional branch on `(~isempty(jmatrix{spin_a,spin_c}))&&(~ismember([spin_labels{1} '_' spin_labels{2} '_' spin_labels{3}],allowed_collision…`.
- Line 306: `for` loop over `n=1:size(subgraphs,1)`.
- Line 307: conditional branch on `any(sum(proxmatrix(subgraphs(n,:),subgraphs(n,:)))==4)`.
- Line 535: `for` loop over `n=1:size(subgraphs,1)`.

### Key state/data transformations

- Lines 47: computes `jmatrix` using `jmatrix=cell(numel(coords),numel(coords))`.
- Lines 50: computes `numbers` using `numbers=1:numel(coords)`.
- Lines 53: computes `proxmatrix` using `proxmatrix=false(numel(coords),numel(coords))`.
- Lines 61: computes `subgraphs` using `subgraphs=dfpt(sparse(proxmatrix),2)`.
- Lines 64: computes `J_NH` using `J_NH=-90.0; J_CN=-15.0; J_CH=140.0; J_CC=35.0`.
- Lines 67: computes `pairs_database` using `pairs_database={`.
- Lines 98: computes `spin_labels` using `spin_labels=pdb_id(subgraphs(n,:))`.
- Lines 99: computes `spin_numbers` using `spin_numbers=numbers(subgraphs(n,:))`.
- Lines 100: computes `spin_resnames` using `spin_resnames=aa_typ(subgraphs(n,:))`.
- Lines 101: computes `spin_resnums` using `spin_resnums=aa_num(subgraphs(n,:))`.
- Lines 104: computes `[spin_labels,index]` using `[spin_labels,index]=sort(spin_labels)`.
- Lines 118: computes `spec_string` using `spec_string=pairs_database(strcmp([spin_labels{1} '_' spin_labels{2}],pairs_database(:,1)),:)`.
- Lines 135: computes `jmatrix{spin_numbers(1),spin_numbers(2)}` using `jmatrix{spin_numbers(1),spin_numbers(2)}=spec_string{2}`.
- Lines 152: computes `J_CCC` using `J_CCC=-1.2`.
- Lines 153: computes `J_CCH` using `J_CCH=0`.
- Lines 154: computes `J_CCN` using `J_CCN=7.0`.
- Lines 155: computes `J_HCH` using `J_HCH=-12.0`.
- Lines 156: computes `J_HNH` using `J_HNH=-1.5`.

### Local helper functions

- Line 633: `grumble()` — `function grumble(aa_num,aa_typ,pdb_id,coords)`.
  - Representative operation: `if ~isnumeric(aa_num)`.
  - Representative operation: `error('aa_num must be a vector of positive integers.')`.

## Parameters / inputs

- aa_num -a vector of amino acid numbers
- aa_typ -a cell array of amino acid types
- pdb_id -a cell array of PDB atom identifiers
- coords -a cell array of coordinate vectors

## Outputs

- jmatrix -a matrix of J-couplings
- Database nomenclature is as follows:
- 1. Atoms in the subgraph descriptor are listed alphabetically to make
- the descriptors unique.
- 2. The four numbers refer to the bonding order, e.g. 1, 3, 2, 4 means
- that the first atom in the descriptor is bonded to the third, which
- is bonded to the second, which is bonded to the fourth. The coupling
- in this case is between atom 1 and atom 4 in the descriptor.
- Note: these J-couplings should be considered approximate. For accurate
- protein work you must supply your own J-couplings.
- Note: this is an auxiliary function that is called by protein.m protein
- import module. Direct calls are discouraged.

## Implementation structure

- Assigns J-couplings from literature values and Karplus curves. Syntax:
- jmatrix=guess_j_pro(aa_num,aa_typ,pdb_id,coords)
- aa_num -a vector of amino acid numbers
- aa_typ -a cell array of amino acid types
- pdb_id -a cell array of PDB atom identifiers
- coords -a cell array of coordinate vectors
- jmatrix -a matrix of J-couplings
- Database nomenclature is as follows:
- 1. Atoms in the subgraph descriptor are listed alphabetically to make
- the descriptors unique.
- 2. The four numbers refer to the bonding order, e.g. 1, 3, 2, 4 means
- that the first atom in the descriptor is bonded to the third, which

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `false()`, `proxmatrix()`, `dfpt()`, `pdb_id()`, `subgraphs()`, `numbers()`, `aa_typ()`, `aa_num()`, `spin_numbers()`, `spin_resnames()`, `spin_resnums()`, `isscalar()`, `num2str()`, `pairs_database()`, `strcmp()`.
