% Test: interfaces/pdb_residues | PDB residue numbers
% Tests protein PDB import across the three- to four-digit residue boundary.
% Syntax:
%
%                     result=test_pdb_residues()
%
% Outputs:
%
%     result - checks for residue and atom identity and coordinates
%
% talos@spindynamics.org

function result=test_pdb_residues()

% Describe the import contract
result=new_test_result('interfaces/pdb_residues','PDB residue numbers',...
                       'adjacent chain and residue fields are read independently.');

% Use consecutive records from deposited 1IRK and one space-delimited record
file_name=[tempname '.pdb'];
cleanup=onCleanup(@()delete(file_name));
lines=["ATOM    147  CD2 LEU A 999      -8.329  74.278  16.744  1.00 17.16";...
       "ATOM    148  N   ARG A1000      -4.096  75.622  19.832  1.00 19.83";...
       "ATOM 149 CA arg A 1000 -3.579 75.879 21.170 1.00 24.60"];
writelines(lines,file_name);

% Check all five outputs in record order
[res,types,names,xyz,serials]=read_pdb_pro(file_name,1);
result=test_true(result,'atom identity',isequal(serials,[147;148;149])&&...
                 isequal(res,[999;1000;1000])&&isequal(types,{'LEU';'ARG';'ARG'})&&...
                 isequal(names,{'CD2';'N';'CA'}),'retain the identities in all three records');
result=test_true(result,'coordinates',isequal(xyz,{[-8.329 74.278 16.744];...
                                                   [-4.096 75.622 19.832];...
                                                   [-3.579 75.879 21.170]}),...
                 'retain the coordinates of each atom');

% Keep both alternate locations distinguishable during BMRB residue matching
lines=["ATOM    547  CG AMET A1051      12.090  68.381   9.108  0.50 36.86";...
       "ATOM    548  CG BMET A1051      12.094  68.360   9.084  0.50 34.17"];
writelines(lines,file_name);
[res,types]=read_pdb_pro(file_name,1);
result=test_true(result,'alternate locations',isequal(res,[1051;1051])&&...
                 isequal(types,{'AMET';'BMET'}),'do not collapse alternate locations onto MET');

end


