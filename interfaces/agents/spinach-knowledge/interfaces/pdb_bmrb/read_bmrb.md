# interfaces/pdb_bmrb/read_bmrb.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/pdb_bmrb/read_bmrb.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=Read_bmrb.m)

`[aa_num,aa_typ,pdb_id,chemsh]=read_bmrb(bmrb_file_name)` reads the Spinach-supported BMRB assignment text format and returns one entry per accepted record. The filename must be a MATLAB character array; there are no optional arguments. The parser replaces tabs with spaces, scans each line using eight space-delimited fields, and accepts a line only when all eight fields parse as nonempty. This is a specific field layout, not a general-purpose reader for arbitrary BMRB formats.

For N accepted assignments, aa_num is an N×1 numeric vector, aa_typ and pdb_id are N×1 cell arrays of strings, and chemsh is an N×1 numeric vector of chemical shifts in ppm. Residue types are converted to uppercase. Values are taken from fields 2, 3, 4 and 6; fields 1, 5, 7 and 8 must nevertheless be present for the record to be retained. An empty or wholly unparseable assignment file raises an error stating that the data format is not canonical.

The source header's output list calls the shift vector chems, but the function signature and implementation use chemsh. Direct calls are discouraged in favor of the protein import workflow described in the protein HOWTO.
