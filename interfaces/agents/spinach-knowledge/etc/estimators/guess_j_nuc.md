# etc/estimators/guess_j_nuc.m

`jmatrix=guess_j_nuc(nuc_num,nuc_typ,pdb_id,coords)`

Estimates RNA one-, two-, and three-bond nuclear-spin J couplings in Hz from nucleotide and atom labels and coordinates. Couplings are selected from the routine's built-in literature and generic-value tables; three-bond entries use Karplus curves. This is a coordinate-based estimator, not a general coupling calculator.

## Inputs and output

- `nuc_num`: nucleotide number for each input atom.
- `nuc_typ`: corresponding nucleotide-type cell array.
- `pdb_id`: corresponding PDB atom-identifier cell array.
- `coords`: corresponding coordinate-vector cell array.
- `jmatrix`: square cell array indexed by input atoms. Assigned elements contain J values in Hz; cells with no supported assignment remain empty.

The input arrays must have matching element counts. The routine infers connectivity using a 1.55 coordinate-distance cutoff, enumerates connected pairs, triples, and four-atom paths, and assigns supported one-, two-, or three-bond couplings from its tables. Three-bond values are evaluated from the stored Karplus coefficients and the path's dihedral angle. The source includes generic or provisional table values as well as literature-derived values, so an estimate should not be mistaken for a measured coupling.

## Table nomenclature and coverage

For a four-atom subgraph descriptor, atom names are listed alphabetically so the descriptor is unique. The four numbers specify the bonded path: for example, `1, 3, 2, 4` means atom 1 is bonded to atom 3, then atom 2, then atom 4; the reported coupling is between atoms 1 and 4. Only patterns present in the built-in pair, triple, and quadruple tables are assigned. Unsupported patterns are not guessed and produce warnings; the routine also prints assignment summaries.

The source validates the input container types and that their element counts match, but does not fully validate coordinate-vector contents. Values and coverage are limited to the hard-coded database and the detected connectivity.

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=guess_j_nuc.m).
