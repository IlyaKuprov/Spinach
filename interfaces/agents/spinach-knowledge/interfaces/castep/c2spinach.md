# interfaces/castep/c2spinach.m

- Signature: `props=c2spinach(file_name)`
- Source: [interfaces/castep/c2spinach.m](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/castep/c2spinach.m)
- Wiki: [c2spinach.m](https://spindynamics.org/wiki/index.php?title=c2spinach.m)

## Purpose and input

Parses a CCP-NC magres v1.0 text file (the source documents CASTEP and other codes) into a MATLAB structure keyed to the atom order in the file. `file_name` must be a character array. The parser removes trailing `#` comments, trims lines, and requires exactly one balanced `[atoms]` block and one `[magres]` block; angle-bracket block tags are also recognised. Other blocks, such as `[magres_old]`, are not used.

Atom rows contain the species, a label, a numeric index, and three coordinates. Atom label/index pairs must be unique. The parser accepts the standard units `Angstrom` for `atom`, `ppm` for `ms`, `au` for `efg`, and `10^19.T^2.J^-1` for `isc` when corresponding units records occur; a different declared unit for these record types is rejected. The code does not perform a coordinate-unit conversion.

## Returned structure

- `props.filename`: the input filename.
- `props.symbols`: a 1-by-`natoms` cell array of atomic symbols in atom-row order.
- `props.std_geom`: an `natoms`-by-3 coordinate matrix in Angstrom; `props.natoms` is the atom count.
- `props.cst`: when at least one `ms` record is present, a 1-by-`natoms` cell array of chemical-shielding tensors in ppm, relative to the bare nucleus. Each tensor is assembled as a 3-by-3 matrix in the printed component order: rows `xx xy xz`, `yx yy yz`, `zx zy zz`.
- `props.efg`: when at least one `efg` record is present, a 1-by-`natoms` cell array of 3-by-3 EFG tensors in atomic units, with the same printed-component ordering.
- `props.k_couplings`: when one or more `isc` records are present, an `natoms`-by-`natoms` symmetric matrix of isotropic reduced couplings in Hz, using the `gparse` convention (magres K multiplied by `mu_N^2/h`).

A missing site tensor is represented by an empty cell; if no record of a tensor type occurs anywhere, its field is omitted. For each `isc` record the code averages the three diagonal tensor components, multiplies by `1e19*(5.0507837461e-27)^2/6.62607015e-34`, and stores that value in both atom-pair directions. Self-couplings are skipped (leaving a zero diagonal); repeated or reversed pair records are averaged. CASTEP's glued label/index form, such as `O100`, is resolved as well as a label and index in separate tokens in magnetic-resonance records; unresolved or ambiguous atom matches are errors.

## Parsing guardrails and consumers

Malformed atom, tensor, or coupling rows; duplicate atom keys; repeated `ms` or `efg` records for one atom; invalid block structure; unrecognised atom matches; and unsupported declared units for these record types raise errors. The argument validator checks the character-array type, while file opening and parsing are performed directly by the routine.

The source documents the reduced-coupling output as input to `g2spinach`, which converts it to isotope-specific J couplings. Existing repository notes also identify `gparse`/`oparse` conventions and consumers in `examples/nmr_solids/case_studies/mathies_14n_13c/*`, `examples/nmr_solids/case_studies/mathies_carbonate/*`, `examples/visualisation/efg_silicate.m` (through `efg_display`), `g2spinach`, and `tests/interfaces/test_c2spinach.m`; those paths are retained as navigation cues, not claims that the examples were run for this note.
