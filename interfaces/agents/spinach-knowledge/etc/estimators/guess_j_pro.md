# etc/estimators/guess_j_pro.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/estimators/guess_j_pro.m`
- Signature: `jmatrix=guess_j_pro(aa_num,aa_typ,pdb_id,coords)`
- Total lines: 653

## Purpose

Assigns J-couplings from literature values and Karplus curves. Syntax: jmatrix=guess_j_pro(aa_num,aa_typ,pdb_id,coords)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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
