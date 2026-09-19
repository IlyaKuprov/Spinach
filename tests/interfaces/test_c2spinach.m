% Tests the magres file parser on the CASTEP files bundled with the
% examples and on small synthetic files. Syntax:
%
%                    result=test_c2spinach()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that tensors are attached to atoms by label and
% index rather than by the order of the records, that the glued
% label and index tokens CASTEP prints for three-digit indices are
% resolved, that isc records become K-couplings in gparse units,
% and that malformed or ambiguous files are refused.
%
% ilya.kuprov@weizmann.ac.il

function result=test_c2spinach()

% Announce the test target
fprintf('TESTING: magres file parser\n');

% State the physical target of the test
result=new_test_result('interfaces/c2spinach',...
                       'magres file parser',...
                       'magnetic parameters are read from magres files and attached to the correct atoms.');

% Locate the example files bundled with Spinach
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));

% Perchlorate, 248 atoms with glued label and index tokens
props=c2spinach(fullfile(root_dir,'examples','standard_systems','castep','perchlorate.magres'));
result=test_true(result,'atom count',props.natoms==248,...
                 'the atom records set the length of every per-atom output');
result=test_true(result,'complete tensor lists',all(~cellfun(@isempty,props.cst))&&all(~cellfun(@isempty,props.efg)),...
                 'every atom, including those with glued records, receives its tensors');
result=test_close(result,'glued record value',props.cst{100}(1,1),26.187706893612472,1e-12,1e-12,...
                  'the record "ms H100" lands on atom H 100');

% Synthetic file with reordered and partial records, a zero-padded index, angle bracket tags, comments, and glued isc endpoints
file_name=[tempname() '.magres'];
file_id=fopen(file_name,'w');
fprintf(file_id,'%s\n','#$magres-abinitio-v1.0','<atoms>','units atom Angstrom',...
                       'atom H H 001 0.0 0.0 0.0 # first','atom C C1 1 1.0 0.0 0.0',...
                       'atom N N 100 0.0 1.0 0.0','</atoms>','<magres>','units ms ppm',...
                       'ms N100 30 0 0 0 30 0 0 0 30','ms H 1 10 1 2 3 10 4 5 6 10',...
                       'units isc 10^19.T^2.J^-1','isc C1 1 N100 3 0 0 0 3 0 0 0 3',...
                       'isc N100 C1 1 6 0 0 0 6 0 0 0 6','isc N100 N100 1 0 0 0 1 0 0 0 1','</magres>');
fclose(file_id); props=c2spinach(file_name); delete(file_name);
result=test_true(result,'record order',isequal(props.cst{3},30*eye(3))&&isempty(props.cst{2}),...
                 'tensors follow the atom label and index, not the order of the records');
result=test_true(result,'zero padded index',isequal(props.cst{1},[10 1 2; 3 10 4; 5 6 10]),...
                 'the printed component order is kept and the index 001 matches the index 1');
k_factor=1e19*(5.0507837461e-27)^2/6.62607015e-34;
result=test_close(result,'averaged K-coupling',props.k_couplings(2,3),4.5*k_factor,1e-9,1e-9,...
                  'the isotropic parts of K(A,B) and K(B,A) are averaged and scaled into gparse units');
result=test_true(result,'no self-coupling',all(diag(props.k_couplings)==0),...
                 'self-coupling records do not reach the diagonal');

% Ambiguous glued token: label C1 index 1 and label C index 11
file_id=fopen(file_name,'w');
fprintf(file_id,'%s\n','[atoms]','atom C C1 1 0 0 0','atom C C 11 1 0 0','[/atoms]',...
                       '[magres]','ms C11 1 0 0 0 1 0 0 0 1','[/magres]');
fclose(file_id);
try
    c2spinach(file_name); refused=false;
catch err
    refused=contains(err.message,'exactly one atom');
end
delete(file_name);
result=test_true(result,'ambiguous glued record',refused,...
                 'a glued token that fits two atoms is an error, not a guess');

end

