# etc/estimators/guess_j_pro.m

- MATLAB implementation: [etc/estimators/guess_j_pro.m](https://github.com/IlyaKuprov/Spinach/blob/main/etc/estimators/guess_j_pro.m)

`jmatrix=guess_j_pro(aa_num,aa_typ,pdb_id,coords)`

Estimates protein J-couplings from atom labels, residue labels, and coordinates. This is an auxiliary routine called by the `protein.m` protein-import module; direct calls are discouraged. The estimates are approximate, and the source explicitly advises supplying your own couplings for accurate protein work.

## Inputs and output

- `aa_num`: numeric vector of residue numbers, one per atom.
- `aa_typ`: cell array of residue/amino-acid type strings, aligned elementwise with the other inputs.
- `pdb_id`: cell array of PDB atom-identifier strings, aligned with the coordinates.
- `coords`: cell array of atom coordinate 3-vectors, one per atom.
- `jmatrix`: preallocated `numel(coords)`-by-`numel(coords)` cell array. Recognised scalar couplings are written at the endpoint atom indices; unassigned cells remain empty. The code does not explicitly mirror an assignment into the reverse-index cell, so callers should not assume symmetry from this routine alone. Reported coupling values are in Hz.

The checker verifies the input classes (numeric for `aa_num`; cells for the other arrays) and equal element counts. Although its error text says positive-integer residue numbers and 3-vectors, it does not actually validate integer/positive values, cell contents, vector shape, atom-label validity, or coordinate units. The geometry code uses a strict distance cutoff of `< 1.60` in the supplied coordinate units; The function does not state or enforce a coordinate unit. Supply mutually aligned atom-level arrays and coordinates in the expected scale.

## How the estimates are made

1. It constructs a proximity matrix by testing every coordinate pair with the Euclidean norm and the 1.60 cutoff. It treats this coordinate-derived graph as connectivity; it does not consume an explicit bond list.
2. Using `dfpt` on that graph, it enumerates connected subgraphs for two-, three-, and four-atom coupling paths. Atom labels are sorted for database lookup. In a descriptor, the index sequence gives the path's bond order: for example `1, 3, 2, 4` means descriptor atom 1—3—2—4, and the coupling endpoints are atoms 1 and 4. The three-atom database similarly specifies the coupling from its first to third atom.
3. One-bond pairs are looked up in an explicit table with backbone-specific values and generic defaults `J_NH=-90.0`, `J_CN=-15.0`, `J_CH=140.0`, and `J_CC=35.0` (Hz). Some backbone entries override these defaults, e.g. CA–CB +34.9, CA–HA +143.5, CA–N −10.7, C–CA +52.5, C–N −14.4, and H–N −93.3 Hz.
4. Three-atom paths use a separate table of scalar values. Its generic constants include `J_CCC=-1.2`, `J_CCH=0`, `J_CCN=-7.0`, `J_HCH=-12.0`, and `J_HNH=-1.5` Hz; several other listed generic terms are zero. The table includes special tryptophan-side-chain values identified in the source as GIAO DFT M06/cc-pVTZ data in PCM water. The source notes two-/three-bond collisions for HIS and PRO and permits the descriptor `H_H_N` as an explicit collision exception; other duplicate assignments raise an error rather than silently overwrite a value.
5. Four-atom paths are filtered to remove T-shaped subgraphs, then assigned through a table of path-specific Karplus coefficients. The evaluated relation is `J = A cos²(theta) + B cos(theta) + C`, where `theta` is the dihedral returned by `dihedral`, and MATLAB's `cosd` makes the angle degrees. Generic coefficient triples (A, B, C; Hz) include aliphatic C–C–C–C [4.48, 0.18, −0.57], aliphatic C–C–C–H [3.83, −0.90, 3.81], and aliphatic H–C–C–H [4.22, −0.50, 4.50]. Aromatic defaults in the source are [0, 0, 5.00], [0, 0, 8.20], and [0, 0, 10.2], respectively. Backbone/side-chain paths and histidine/proline rings also have specific table entries; the latter are identified as GIAO DFT M06/cc-pVTZ data.

The function prints a one-bond summary and estimated assignments. Isolated spins and unrecognised pair/triple/quadruple patterns produce warnings; duplicate table matches or conflicting assignments can error. These values are a table-and-geometry approximation, not a general quantum-chemical calculation or guarantee of complete coverage.

## Source references

- [Function wiki page](https://spindynamics.org/wiki/index.php?title=guess_j_pro.m)
- The source associates the generic aliphatic Karplus sets with [Carbohydrate Research DOI 10.1016/j.carres.2007.02.023](https://doi.org/10.1016/j.carres.2007.02.023), [Angewandte Chemie DOI 10.1002/anie.198204491](https://doi.org/10.1002/anie.198204491), and [JACS DOI 10.1021/ja00901a059](https://doi.org/10.1021/ja00901a059). It cites Vuister, “J-couplings. Measurement and Usage in Structure Determination,” for chi1 side-chain entries, and DOI 10.1021/ja00111a021 for backbone psi entries.
- The source labels aromatic defaults only as “[literature reference]” and supplies no fuller citation there. It gives method labels, but no bibliographic citation, for the GIAO DFT special-case data; none is inferred here.
