# etc/estimators/guess_j_pro.m

`jmatrix=guess_j_pro(aa_num,aa_typ,pdb_id,coords)`

Estimates approximate protein J-couplings from literature values and Karplus curves. This auxiliary function is called by the `protein.m` protein-import module; direct calls are discouraged.

## Inputs and output

- `aa_num`: vector of amino-acid numbers.
- `aa_typ`: cell array of amino-acid types.
- `pdb_id`: cell array of PDB atom identifiers.
- `coords`: cell array of coordinate vectors.
- `jmatrix`: square cell array of J-couplings in Hz; unassigned entries remain empty.

The four inputs must have the same length.

## Method and limitations

The function infers bonds from interatomic distances below 1.60 Å, then examines connected atom pairs, triples, and quadruples. One- and two-bond couplings are selected by matching atom labels to tabulated values. Three-bond couplings are estimated with `A*cosd(theta)^2+B*cosd(theta)+C`, using a coordinate-derived dihedral. T-shaped quadruples are excluded. Isolated atoms and unrecognized groups produce warnings; conflicting assignments can produce errors.

Database descriptors list atom labels alphabetically. Their four indices specify bonding order: `1, 3, 2, 4` means descriptor atom 1 bonds to atom 3, then atom 2, then atom 4; the coupling is between atoms 1 and 4.

These estimates are approximate. For accurate protein work, supply your own J-couplings.