# interfaces/pdb_bmrb/read_pdb_pro.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/pdb_bmrb/read_pdb_pro.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=read_pdb_pro.m)

`[aa_num,aa_typ,pdb_id,coords,pdb_ser]=read_pdb_pro(pdb_file_name,mod_id)` reads protein ATOM records from a selected PDB model. Both arguments are required: the filename must be a MATLAB character array, and mod_id must be a finite, real, positive integer scalar. A model-free file is read when mod_id=1 (the parser reports that it is reading a single-model file); there is no default for the second argument.

The reader scans for a matching MODEL record, then parses ATOM lines until ENDMDL or end of file. Standard records use fixed PDB columns, including residue sequence number in columns 23–26 and Cartesian coordinates in columns 31–54. Historical space-delimited records with all ten expected fields are also accepted. For N accepted atoms, aa_num and pdb_ser are N×1 numeric vectors; aa_typ and pdb_id are N×1 cell arrays of strings; and coords is an N×1 cell array of three-element coordinate vectors in ångströms. Alternate-location identifiers prefix the residue label (for example, AMET and BMET), so resolve alternate locations before matching ordinary BMRB residue names. Residue types are uppercased; examples are 'TYR' and atom 'HE2'. The serial-number output is the PDB atom serial from the input record.

If a requested model is absent, the function has no dedicated not-found error; it can return no parsed records.
