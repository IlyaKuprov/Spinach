# interfaces/pdb_bmrb/read_pdb_nuc.m

[MATLAB source](https://github.com/IlyaKuprov/Spinach/blob/main/interfaces/pdb_bmrb/read_pdb_nuc.m) · [Spin Dynamics Wiki](https://spindynamics.org/wiki/index.php?title=Read_pdb_nuc.m)

`[res_num,res_typ,pdb_id,coords]=read_pdb_nuc(pdb_file_name)` reads nucleic-acid atom records from a PDB text file. The filename must be a MATLAB character array; there are no optional arguments. It scans lines using the space-delimited ATOM record pattern and retains a line only when all ten parsed fields are nonempty. Thus it does not promise to accept every valid PDB record variant or HETATM records.

For N accepted records, res_num is an N×1 numeric vector, res_typ and pdb_id are N×1 cell arrays of strings, and coords is an N×1 cell array whose entries are three-element coordinate vectors in ångströms. Residue types are converted to uppercase; examples include residue 'GUA' and atom 'C1P'. The parser reads the chain identifier from character 22 when present, but does not return it; multiple distinct chain identifiers among accepted records cause an error. It does not select or validate a single model, so prepare input as a single-model, single-chain file when that is required by the downstream workflow.
