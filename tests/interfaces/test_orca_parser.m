% Tests the ORCA log parser on the logs bundled with the examples. Syntax:
%
%                    result=test_orca_parser()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks version detection, the mapping of ORCA's zero-based
% and incompletely printed nucleus numbering onto the atoms of the
% coordinate table, and the numerical values of the g-tensor and of
% the hyperfine tensors against the quantities that ORCA itself
% prints separately in the same logs.
%
% ilya.kuprov@weizmann.ac.il

function result=test_orca_parser()

% Announce the test target
fprintf('TESTING: ORCA log parser\n');

% State the physical target of the test
result=new_test_result('interfaces/orca_parser',...
                       'ORCA log parser',...
                       'magnetic parameters are read from ORCA logs and attached to the correct atoms.');

% Locate the example logs bundled with Spinach
root_dir=fileparts(fileparts(fileparts(mfilename('fullpath'))));

% Methyl radical, a vacuum DFT calculation with a g-tensor and hyperfines
props=oparse(fullfile(root_dir,'examples','esr_liq_pulsed',...
                      'data_import','orca_methyl_radical.out'));
result=test_true(result,'version detection',strcmp(props.orca_version,'4.1.1'),...
                 'the parser branch is chosen from the ORCA version banner');
result=test_close(result,'methyl radical g-tensor',props.g_tensor.matrix(1,1),...
                  2.0022269,1e-7,1e-7,...
                  'the g-matrix is read from the electronic g-matrix section');
proton_hfc=gauss2mhz(props.hfc.full.matrix{2});
result=test_close(result,'methyl radical hyperfine',trace(proton_hfc)/3,...
                  -65.1341,1e-4,1e-4,...
                  'the isotropic hyperfine coupling matches the A(iso) printed by ORCA');

% Copper porphyrin, hyperfine tensors printed for the protons only
props=oparse(fullfile(root_dir,'examples','visualisation','porphyrine.out'));
tensor_atoms=find(~cellfun(@isempty,props.hfc.full.matrix));
proton_atoms=find(strcmp(props.symbols,'H'));
result=test_true(result,'skipped nucleus numbering',isequal(tensor_atoms,proton_atoms),...
                 'ORCA numbers printed nuclei from zero and skips the ones it was not asked about');
result=test_true(result,'atom count',props.natoms==37,...
                 'the coordinate table sets the length of every per-atom output');

end

