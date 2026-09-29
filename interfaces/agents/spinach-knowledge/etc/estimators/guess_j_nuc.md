# etc/estimators/guess_j_nuc.m

- MATLAB implementation: [etc/estimators/guess_j_nuc.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/estimators/guess_j_nuc.m)

`jmatrix=guess_j_nuc(nuc_num,nuc_typ,pdb_id,coords)`

Estimate RNA nuclear-spin J couplings in Hz from atom labels, nucleotide labels, and coordinates. The result is a square cell matrix indexed by input-atom positions; supported assignments occupy the table-designated coupling endpoints, and the routine does not mirror each value into the opposite matrix element. Empty cells mean no assignment was made, not a zero coupling.

## Inputs

Each input describes the same atom list: `nuc_num` is the nucleotide number, `nuc_typ` the nucleotide type, `pdb_id` the PDB atom identifier, and `coords` the atom coordinate vector. The labels and coordinates are cell arrays and all four inputs must have equal element counts. The source documents nucleotide numbers as positive integers and coordinates as 3-vectors, but its validator checks only that `nuc_num` is numeric, that the other inputs are cells, and that their element counts match; it does not validate integer values, coordinate shape, or coordinate units.

## How assignments are made

Connectivity is inferred by joining every pair of coordinates whose Euclidean separation is less than 1.55 in the supplied coordinate units. The routine enumerates connected pairs, triples, and four-atom paths and matches their atom labels against RNA-specific and generic tables. Consequently, the coordinate scale and atom naming must be compatible with the source's cutoff and nomenclature; this is not a general J-coupling calculator.

- **One bond:** a pair-label lookup uses literature entries where present and generic values otherwise. Generic examples in the source are `J_NH = −86`, `J_CN = −20`, `J_CH = 180`, and `J_CC = 65 Hz`.
- **Two bonds:** a three-atom label key selects a coupling and specifies which ordered connected path the key represents. Examples of generic values are `J_CCC = −1.2`, `J_CCN = −7`, `J_HCH = −12`, and `J_HNH = 10 Hz`; several other generic combinations are set to zero. The source marks some ring values as provisional or to be replaced by DFT values.
- **Three bonds:** a four-atom key specifies the bond path and coupling endpoints, together with Karplus coefficients. For the path dihedral `theta`, the stored parameters are evaluated as `A*cosd(theta)^2 + B*cosd(theta) + C`. The table is selective and many entries have zero coefficients, so a Karplus calculation is made only for a matched listed path.

Atom labels in each database key are alphabetised to make keys unique; the associated index sequence restores the bonding order and identifies the coupling endpoints. Unsupported connected pairs or paths are warned about and left unassigned; duplicate database records are errors, and conflicting assignments can also raise an error. Four-atom branched (T-shaped) subgraphs are excluded before the path lookup.

The table mixes literature-specific values with generic/provisional values, so results are estimates conditioned on the supplied topology and hard-coded coverage—not experimental measurements or a complete RNA coupling model.

[Spin Dynamics Wiki page](https://spindynamics.org/wiki/index.php?title=Guess_j_nuc.m).
