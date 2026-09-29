# interfaces/pdb_bmrb/read_pdb_pro.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/pdb_bmrb/read_pdb_pro.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=read_pdb_pro.m)

`[aa_num,aa_typ,pdb_id,coords,pdb_ser]=read_pdb_pro(pdb_file_name,mod_id)` reads protein ATOM records from a selected PDB model. Both arguments are required: the filename must be a MATLAB character array, and mod_id must be a finite, real, positive integer scalar. A model-free file is read when mod_id=1 (the parser reports that it is reading a single-model file); there is no default for the second argument.

The reader scans for a matching MODEL record, then parses space-delimited ATOM lines until ENDMDL or end of file. An accepted line must provide all ten expected fields. For N accepted atoms, aa_num and pdb_ser are N×1 numeric vectors; aa_typ and pdb_id are N×1 cell arrays of strings; and coords is an N×1 cell array of three-element coordinate vectors in ångströms. Residue types are uppercased; examples are 'TYR' and atom 'HE2'. The serial-number output is the PDB atom serial from the input record.

If a requested model is absent, the function has no dedicated not-found error; it can return no parsed records. The supported input is the implementation's expected space-delimited ATOM layout, rather than every possible PDB representation.
