# etc/estimators/guess_j_nuc.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/etc/estimators/guess_j_nuc.m`
- Signature: `jmatrix=guess_j_nuc(nuc_num,nuc_typ,pdb_id,coords)`
- Total lines: 465

## Purpose

RNA assignments of J-couplings from literature values and Karplus cur- ves. Syntax: jmatrix=guess_j_nuc(nuc_num,nuc_typ,pdb_id,coords)

## Physical / mathematical content

- This file belongs to the `etc` part of Spinach. Its role should be read together with nearby files in the same directory, which usually share a common physical regime or infrastructure purpose.

## Numerical / algorithmic content

- The file contains an explicit `grumble(...)` validator, which is Spinach convention for front-loading dimension, type, and regime checks before expensive linear-algebra work begins.
- The file also defines local helper function(s): `numel()`. This usually means the public entry point is supported by tightly coupled validation or helper logic kept private to the file.

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
