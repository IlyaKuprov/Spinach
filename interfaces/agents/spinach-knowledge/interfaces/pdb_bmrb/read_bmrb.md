# interfaces/pdb_bmrb/read_bmrb.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/pdb_bmrb/read_bmrb.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=Read_bmrb.m)

`[aa_num,aa_typ,pdb_id,chemsh]=read_bmrb(bmrb_file_name)` extracts chemical shift assignments from NMR-STAR 2.1 and 3.x BMRB files, or from the historical headerless eight-column tables. The filename is a non-empty MATLAB character row vector; there are no optional arguments. Prefer `protein` for the combined PDB, BMRB, and CSA import workflow described in the protein HOWTO.

The tagged reader identifies residue numbers, residue labels, atom names, and shifts by their column tags, rather than assuming positions. Columns may be reordered or optional fields omitted. Comments, quoted values, semicolon-delimited multiline text, and rows wrapped across lines do not change the assignments. Unrelated loops and metadata are not mistaken for shifts. Shift loops may end with `stop_`, another structural token, or the end of the file.

For N accepted assignments, `aa_num` and `chemsh` are N×1 numeric vectors, while `aa_typ` and `pdb_id` are N×1 cell arrays of character vectors. Shifts are in ppm, residue labels are uppercased, and atom labels are retained without renaming. Numeric author residue numbers take precedence over sequence and component-index numbering, with per-row fallback where author numbering is unavailable. Null assignments are omitted; no missing shift is invented. Incomplete rows, missing or duplicated assignment tags, and invalid non-null numerical assignments raise errors. A valid STAR file with no shift assignments returns empty columns.

All entities and chemical shift lists are returned in file order. This four-output interface does not select a chain, resolve repeated assignments across lists, or represent insertion-code author identifiers; non-numeric author numbering falls back to the numeric sequence identifiers. Successful parsing therefore does not establish that a PDB and BMRB file describe the same sequence or numbering. Those consistency checks remain in `protein`.

The legacy headerless path retains the original layout: fields 2, 3, 4, and 6 supply the outputs, and all eight fields must parse. It is a compatibility path, not an alternative interpretation of tagged STAR files.

Format reference: [BMRB assigned chemical shift file help](https://bmrb.io/deposit/shifts_example_help.shtml).
