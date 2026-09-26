# etc/estimators/guess_j_nuc.m

- Signature: `jmatrix=guess_j_nuc(nuc_num,nuc_typ,pdb_id,coords)`

## Purpose

RNA assignments of J-couplings from literature values and Karplus cur- ves. Syntax: jmatrix=guess_j_nuc(nuc_num,nuc_typ,pdb_id,coords)

## Physical / mathematical content

## Numerical / algorithmic content

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
