# interfaces/pdb_bmrb/nuclacid.m

Source: [interfaces/pdb_bmrb/nuclacid.m](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/pdb_bmrb/nuclacid.m)
Wiki: [Nuclacid.m](https://spindynamics.org/wiki/index.php?title=Nuclacid.m)

## Interface and inputs

`[sys, inter] = nuclacid(pdb_file, shift_file, options)` imports a nucleic-acid structure and chemical shifts for Spinach. Both file names must be MATLAB character arrays. The required `options.noshift` character value is `'keep'` or `'delete'`; `options.deut_list` must be a cell array. The shift input is ASCII with tab-delimited residue number, atom ID, and shift fields (the parser reads each line as numeric, text, numeric). See `example.txt` in `examples/nmr_nucleic`.

The PDB reader supplies residue numbers/types, atom IDs, and coordinates. Before matching, the importer removes the named atom IDs `O5'`, `O4'`, `O4`, `O6`, `O2`, `O2'`, `O3'`, `O1P`, `O2P`, `P`, `H2'`, `H5T`, and `HO'2`; these atoms do not appear in the returned simulation. Shift rows are selected by residue number and matched on exact atom ID; the match does not additionally compare residue type. When duplicate rows match a residue number and atom ID, the first selected shift is used. PDB atom-name primes are changed to `p` before J-coupling estimation (for example, an apostrophe in an atom ID becomes `p`).

Only retained atom IDs beginning with `H`, `C`, `N`, or `P` are mapped, respectively, to `1H`, `13C`, `15N`, or `31P`; another initial atom character raises an error. Couplings are estimated by `guess_j_nuc` from the retained residue numbers, residue types, normalised atom IDs, and coordinates. The importer does not estimate CSA here.

## Missing shifts and deuteration

For an unassigned retained atom, `options.noshift = 'keep'` assigns values from `linspace(-1,0,nmissing)` ppm to the missing shifts; `'delete'` removes those atoms from both the spin list and both dimensions of the coupling array. The example deuteration key is `'ADE:H2pp'`. Keys are matched as residue type, colon, and normalised atom ID (not the numeric residue index). A requested deuteration must identify a proton mapped as `1H`; otherwise the function errors. For accepted entries, the isotope becomes `2H` and the corresponding coupling row and column are scaled by `spin('2H')/spin('1H')`.

## Outputs

The output order follows the retained PDB atom order after exclusions and, for `'delete'`, after removal of unassigned atoms.

- `sys.isotopes`: row cell array of isotope labels, one per spin.
- `sys.labels`: row cell array of labels in the form `RESNAME(residue_number):atom_id`.
- `inter.coordinates`: Nspins x 3 numeric coordinates in Angstrom.
- `inter.zeeman.scalar`: row cell array of isotropic chemical shifts in ppm.
- `inter.coupling.scalar`: Nspins x Nspins cell array of scalar couplings in Hz.

The returned cell orientations above follow the implementation's row-vector construction, including `inter.zeeman.scalar`; coordinates remain a numeric matrix.
