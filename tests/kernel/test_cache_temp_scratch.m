% Tests cache management in a temporary scratch directory. Syntax:
%
%                    result=test_cache_temp_scratch()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks cacheman() and wipe_cache() against a temporary scratch
% directory, loads shipped small cache tables without generating new kernel
% cache files, and smoke-tests the read-only sniff() integrity pass.
%
% ilya.kuprov@weizmann.ac.il

function result=test_cache_temp_scratch()

% State the cache-management target of the test
result=new_test_result('kernel/cache_temp_scratch',...
                       'Temporary scratch cache management',...
                       'cacheman() and wipe_cache() must only affect Spinach cache records in the configured scratch directory.');

% Create an isolated scratch directory and arrange cleanup
scratch=tempname(tempdir);
mkdir(scratch);
cleanup_obj=onCleanup(@()local_remove_dir(scratch));

% Build a minimal Spinach object pointing at the temporary scratch directory
spin_system.sys.output='hush';
spin_system.sys.scratch=scratch;
spin_system.sys.enable={};
spin_system.sys.disable={};
spin_system.tols.cache_mem=365;

% Ensure cacheman() uses a small process pool instead of auto-starting a large one
current_pool=gcp('nocreate');
if isempty(current_pool)
    parpool('Processes',1);
end

% Create Spinach and non-Spinach scratch files
payload=1;
spinach_file=fullfile(scratch,'spinach_fresh.mat');
ordinary_file=fullfile(scratch,'ordinary_file.mat');
save(spinach_file,'payload');
save(ordinary_file,'payload');

% Check that a fresh Spinach cache file is kept by a long retention horizon
cacheman(spin_system);
result=test_true(result,'cacheman keeps fresh cache file',exist(spinach_file,'file')==2,...
                 'fresh Spinach cache files must survive a long cache_mem horizon');
result=test_true(result,'cacheman ignores ordinary file',exist(ordinary_file,'file')==2,...
                 'files not matching the spinach_* cache pattern must be ignored');

% Check that wipe_cache() deletes Spinach cache files but not unrelated files
wipe_cache(spin_system);
result=test_true(result,'wipe_cache removes Spinach file',exist(spinach_file,'file')==0,...
                 'wipe_cache() must remove matching Spinach cache files in scratch');
result=test_true(result,'wipe_cache keeps ordinary file',exist(ordinary_file,'file')==2,...
                 'wipe_cache() must not remove unrelated scratch files');

% Check that cacheman() removes matching cache directories under a zero horizon
spinach_dir=fullfile(scratch,'spinach_old_dir');
mkdir(spinach_dir);
save(fullfile(spinach_dir,'payload.mat'),'payload');
spin_system.tols.cache_mem=0;
cacheman(spin_system);
result=test_true(result,'cacheman removes Spinach directory',exist(spinach_dir,'dir')==0,...
                 'cacheman() must remove matching out-of-date Spinach cache directories');

% Locate the shipped kernel cache directory
cache_root=fileparts(which('cacheman'));

% Verify that small shipped cache records exist before loading them
result=test_true(result,'st_product_table shipped cache',exist(fullfile(cache_root,'st_product_table_2.mat'),'file')==2,...
                 'the shipped two-level ST product table cache should be available');
result=test_true(result,'bos_product_table shipped cache',exist(fullfile(cache_root,'bos_product_table_2.mat'),'file')==2,...
                 'the shipped two-level bosonic product table cache should be available');
result=test_true(result,'sle_operators shipped cache',exist(fullfile(cache_root,'sle_operators_rank_1.mat'),'file')==2,...
                 'the shipped rank-one SLE operator cache should be available');

% Load small cache-table records and check their dimensions
[st_left,st_right]=st_product_table(2);
result=test_true(result,'st_product_table shape',isequal(size(st_left),[4 4 4])&&isequal(size(st_right),[4 4 4]),...
                 'two-level single-transition product tables must have 4x4x4 dimensions');
[ist_left,ist_right]=ist_product_table(2);
result=test_true(result,'ist_product_table shape',isequal(size(ist_left),[4 4 4])&&isequal(size(ist_right),[4 4 4]),...
                 'spin-half IST product tables must have 4x4x4 dimensions');
[bos_left,bos_right]=bos_product_table(2);
result=test_true(result,'bos_product_table shape',isequal(size(bos_left),[4 4 4])&&isequal(size(bos_right),[4 4 4]),...
                 'two-level bosonic product tables must have 4x4x4 dimensions');
[Lx,Ly,Lz,D,space_basis]=sle_operators(1);
result=test_true(result,'sle_operators dimensions',isequal(size(Lx),[10 10])&&isequal(size(Ly),[10 10])&&...
                 isequal(size(Lz),[10 10])&&isequal(size(D),[5 5])&&isequal(size(space_basis),[10 3]),...
                 'rank-one SLE operators must match the ten-function Wigner basis');

% Smoke-test the read-only integrity sniffer without requiring a pristine tree
old_dir=pwd;
cd_cleanup=onCleanup(@()cd(old_dir));
cd(fileparts(which('sniff')));
sniff_text=evalc('sniff(''none'');');
result=test_true(result,'sniff read-only smoke',~isempty(sniff_text),...
                 'sniff(''none'') should complete and print either an all-clear or fishy-file report');

end


function local_remove_dir(path_name)

% Remove a temporary directory if it still exists
if exist(path_name,'dir')
    rmdir(path_name,'s');
end

end


