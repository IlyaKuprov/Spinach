# interfaces/castep/c2spinach.m

- Source: `/home/kuprov/.openclaw/workspace/Spinach/interfaces/castep/c2spinach.m`
- Signature: `props=c2spinach(file_name)`
- Total lines: 197

## Purpose

Parser for .magres files written by CASTEP and other codes in the CCP-NC magres v1.0 format. Reads the [atoms] and [magres] blocks and returns the geometry and the magnetic resonance tensors, keyed to the atoms in the order in which the atom records appear in the file. Syntax: props=c2spinach(file_name)

## Physical / mathematical content

- CASTEP interface. Recovers the periodic DFT (GIPAW) shielding, electric field gradient, and reduced spin-spin coupling tensors and converts them to the conventions used by `gparse` and `oparse`, so that `g2spinach` consumes the output the same way.
- Shielding tensors are returned in the printed component order (xx xy xz yx yy yz zx zy zz as the rows of the 3x3 matrix), which is the orientation in which Spinach contracts a Zeeman tensor with the spin operator on the left and the field on the right; the antisymmetric part of the shielding therefore enters with the sign the magres standard defines.
- Reduced couplings K (magres units 10^19 T^2 J^-1) are multiplied by mu_N^2/h, which is the unit `gparse` reports Gaussian K-couplings in; `g2spinach` then converts them into J-couplings for the isotopes it is given.

## Numerical / algorithmic content

- Every record is matched with a regular expression that spells out the field count and the number syntax; a record with the right tag and the wrong shape is an error rather than a silently shifted tensor.
- Block tags (`[atoms]`, `[magres]`, also the angle-bracket form) are located first; nested, unbalanced, missing, or repeated [atoms] and [magres] blocks are errors. Everything outside those two blocks, including the `[magres_old]` block CASTEP appends, is ignored.
- Units records for atom, ms, efg, and isc are checked against the standard units and anything else is an error; an absent units record means the standard unit, as the specification prescribes.
- The file contains an explicit `grumble(...)` validator for the input argument and a local key resolver shared by the ms, efg, and isc records.

## Parameters / inputs

- file_name - the name of the *.magres file, a character string

## Outputs

- props.filename - log file name
- props.symbols - atomic symbols, 1 x natoms cell
- props.std_geom - atomic coordinates, natoms x 3, Angstrom
- props.natoms - number of atoms
- props.cst - chemical shielding tensors relative to the bare nucleus in vacuum, ppm, 1 x natoms cell, printed component order
- props.efg - EFG tensors, a.u., 1 x natoms cell
- props.k_couplings - isotropic reduced spin-spin couplings, natoms x natoms, in the units used by gparse.m (magres K times mu_N^2/h, Hz)
- Only the tensors that the file contains are returned; test for their presence with isfield. An atom without a tensor gets an empty cell.

## Internal Spinach / MATLAB structure cues

- Called routines detected from the main body: `grumble()`, `key_atoms()`, `fopen()`, `textscan()`, `fclose()`, `strtrim()`, `regexprep()`, `regexp()`, `cellfun()`, `strcmp()`, `strcat()`, `strsplit()`, `str2double()`, `reshape()`, `unique()`, `isstrprop()`, `ischar()`.
- Consumers in the repository: `examples/nmr_solids/case_studies/mathies_14n_13c/*`, `examples/nmr_solids/case_studies/mathies_carbonate/*`, `examples/visualisation/efg_silicate.m` (through `efg_display`), `g2spinach` for `cst` and `k_couplings`, and `tests/interfaces/test_c2spinach.m`.
