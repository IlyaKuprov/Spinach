% Tests difficult dynamic coverage for includes, integrity, and MEX helpers. Syntax:
%
%              result=test_dynamic_integrity_includes_mex_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test uses direct include execution, read-only integrity probes, and
% temporary-directory fixtures for mutating integrity and MEX helpers. It
% avoids touching production tests, production code, shipped build outputs,
% and repository state.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_integrity_includes_mex_suite()

% State the integrity/include/MEX target of the test
result=new_test_result('kernel/dynamic_integrity_includes_mex',...
                       'Dynamic integrity, include, and MEX helper coverage',...
                       'low-feasibility include scripts and integrity utilities must be exercised without mutating the Spinach tree.');

% Locate canonical Spinach subtrees
spinach_root='/home/kuprov/.openclaw/workspace/Spinach';
includes_dir=fullfile(spinach_root,'kernel','includes');
integrity_dir=fullfile(spinach_root,'kernel','integrity');
mex_dir=fullfile(spinach_root,'kernel','mex');

% Exercise host overrides and GPU guard includes
result=local_test_autoexec(result,fullfile(includes_dir,'autoexec.m'));
result=local_test_gpu_guard(result,fullfile(includes_dir,'start_disallow_gpu.m'),...
                            fullfile(includes_dir,'end_disallow_gpu.m'));

% Start a temporary pool for the pool-dependent include paths
pool_created=local_start_pool_if_needed();
pool_cleanup=onCleanup(@()local_delete_pool(pool_created));

% Exercise parallel profiler and Redfield include paths
result=local_test_parallel_profiler(result,fullfile(includes_dir,'parallel_profiler_start.m'),...
                                    fullfile(includes_dir,'parallel_profiler_report.m'));
result=local_test_redfield_serial(result,fullfile(includes_dir,'redfield_integral_serial.m'));
result=local_test_redfield_async(result,fullfile(includes_dir,'redfield_integral_async.m'));
result=local_test_direct_include_dispatch(result);

% Exercise integrity utilities using read-only and temporary fixtures
result=local_test_existentials(result,fullfile(integrity_dir,'existentials.m'));
result=local_test_exorcise_patrol(result,fullfile(integrity_dir,'exorcise.m'),...
                                  fullfile(integrity_dir,'patrol.m'));
result=local_test_rearm_sniff_fixture(result,fullfile(integrity_dir,'rearm.m'),...
                                      fullfile(integrity_dir,'sniff.m'));
result=local_test_real_sniff(result,integrity_dir);
result=local_test_smack_static(result,fullfile(integrity_dir,'smack.m'));

% Exercise MEX compiler helper without building in the real tree
result=local_test_compile_mex(result,fullfile(mex_dir,'compile_mex.m'));

% Release the temporary pool before returning from the test
clear('pool_cleanup');

end


function result=local_test_autoexec(result,autoexec_file)

% Preserve environment and graphics defaults touched by autoexec
old_computer=getenv('COMPUTERNAME');
old_position=get(groot,'defaultFigurePosition');
old_style=get(groot,'defaultFigureWindowStyle');
old_menu=get(groot,'defaultFigureMenuBar');
old_toolbar=get(groot,'defaultFigureToolbar');
cleanup_obj=onCleanup(@()local_restore_autoexec(old_computer,old_position,...
                                                old_style,old_menu,old_toolbar));

% Check the ELMINSTER process-count override
sys=struct();
setenv('COMPUTERNAME','ELMINSTER');
run(autoexec_file);
result=test_true(result,'autoexec ELMINSTER override',...
                 isequal(sys.parallel,{'processes',128}),...
                 'ELMINSTER should select the documented 128-process profile when the user did not set sys.parallel');

% Check the ALAUNDO GPU-aware override
sys=struct();
sys.enable={'gpu'};
setenv('COMPUTERNAME','ALAUNDO');
run(autoexec_file);
result=test_true(result,'autoexec ALAUNDO GPU override',...
                 isequal(sys.parallel,{'processes',32}),...
                 'ALAUNDO should use four workers per GPU when GPU arithmetic was requested');

% Check the TALOS non-GPU override
sys=struct();
sys.enable={};
setenv('COMPUTERNAME','TALOS');
run(autoexec_file);
result=test_true(result,'autoexec TALOS CPU override',...
                 isequal(sys.parallel,{'processes',56}),...
                 'TALOS should use the documented CPU process profile without GPU arithmetic');

% Check that user parallel settings are not overwritten
sys=struct();
sys.parallel={'processes',4};
sys.enable={'gpu'};
setenv('COMPUTERNAME','TALOS');
run(autoexec_file);
result=test_true(result,'autoexec preserves user parallel setting',...
                 isequal(sys.parallel,{'processes',4}),...
                 'autoexec must not override a user-supplied sys.parallel field');

% Check that unknown hosts leave sys.parallel unset
sys=struct();
setenv('COMPUTERNAME','UNKNOWN_SPINACH_HOST');
run(autoexec_file);
result=test_true(result,'autoexec unknown host no-op',~isfield(sys,'parallel'),...
                 'unknown host names should not invent a parallelisation policy');

end


function result=local_test_gpu_guard(result,start_file,end_file)

% Check GPU removal and restoration around the include pair
spin_system=local_quiet_spin_system();
spin_system.sys.enable={'gpu','mex'};
run(start_file);
result=test_true(result,'start_disallow_gpu removes gpu',...
                 user_wanted_gpu&&(~ismember('gpu',spin_system.sys.enable))&&...
                 ismember('mex',spin_system.sys.enable),...
                 'start_disallow_gpu must remember the request, remove only gpu, and keep unrelated flags');
run(end_file);
result=test_true(result,'end_disallow_gpu restores gpu',...
                 ismember('gpu',spin_system.sys.enable),...
                 'end_disallow_gpu must restore GPU arithmetic when the user had requested it');

% Check the no-GPU branch leaves the enable list unchanged
clear('user_wanted_gpu');
spin_system=local_quiet_spin_system();
spin_system.sys.enable={'mex'};
run(start_file);
run(end_file);
result=test_true(result,'gpu guard no-GPU branch',...
                 (~ismember('gpu',spin_system.sys.enable))&&ismember('mex',spin_system.sys.enable),...
                 'the guard pair must not add gpu when it was not present initially');

% Check that the end include refuses to run without the start include
clear('user_wanted_gpu');
err_text=local_error_text(@()run(end_file));
result=test_true(result,'end_disallow_gpu requires start',...
                 contains(err_text,'must be preceded by start_disallow_gpu'),...
                 'end_disallow_gpu must protect callers from restoring an undefined GPU policy');

end


function result=local_test_parallel_profiler(result,start_file,report_file)

% Build a quiet spin system without the detailed dafuq profiler
spin_system=local_quiet_spin_system();
spin_system.sys.enable={};
spin_system.sys.scratch=tempdir;

% Exercise brief profiler start and report includes
run(start_file);
pause(0.01);
run(report_file);
result=test_true(result,'parallel profiler byte counters',...
                 exist('nbytes','var')&&isnumeric(nbytes)&&numel(nbytes)==2&&...
                 all(isfinite(nbytes)),...
                 'parallel profiler report must compute a finite two-element worker byte counter');
result=test_true(result,'parallel profiler wall time',...
                 exist('walltime','var')&&isnumeric(walltime)&&isscalar(walltime)&&...
                 isfinite(walltime)&&(walltime>=0),...
                 'parallel profiler report must compute a finite non-negative wall time');

% Inspect the detailed internal-profiler branch without executing it
start_src=fileread(start_file);
report_src=fileread(report_file);
result=test_true(result,'parallel profiler dafuq source guard',...
                 contains(start_src,'parallel.internal.profiling.PoolProfiler')&&...
                 contains(report_src,'parProfiler.drainLog()'),...
                 'the detailed dafuq profiler path is internal MATLAB API code and is source-guarded here');

end


function result=local_test_redfield_serial(result,serial_file)

% Build a one-dimensional Redfield fixture with an analytical integral
[spin_system,Q,L0,R,expected]=local_redfield_fixture(); %#ok<ASGLU>
run(serial_file);
result=test_close(result,'serial Redfield include integral',R,expected,1e-10,1e-10,...
                  'the serial include must add the one-dimensional analytical Redfield integral');
result=test_true(result,'serial Redfield clears large inputs',...
                 (~exist('Q','var'))&&(~exist('L0','var')),...
                 'the serial include should clear the large input cells after integration');

end


function result=local_test_redfield_async(result,async_file)

% Build a one-dimensional Redfield fixture with an analytical integral
[spin_system,Q,L0,R,expected]=local_redfield_fixture(); %#ok<ASGLU>
run(async_file);
result=test_close(result,'async Redfield include integral',R,expected,1e-10,1e-10,...
                  'the asynchronous include must reproduce the analytical Redfield integral through parfeval and ValueStore');
result=test_true(result,'async Redfield clears futures and store',...
                 (~exist('F','var'))&&(~exist('store','var')),...
                 'the asynchronous include should clear future and ValueStore handles after assembly');

end


function result=local_test_direct_include_dispatch(result)

% Preserve host and graphics defaults touched by direct include calls
old_computer=getenv('COMPUTERNAME');
old_position=get(groot,'defaultFigurePosition');
old_style=get(groot,'defaultFigureWindowStyle');
old_menu=get(groot,'defaultFigureMenuBar');
old_toolbar=get(groot,'defaultFigureToolbar');
cleanup_obj=onCleanup(@()local_restore_autoexec(old_computer,old_position,...
                                                old_style,old_menu,old_toolbar));

% Execute autoexec by script name in a harmless unknown-host case
sys=struct();
setenv('COMPUTERNAME','UNKNOWN_SPINACH_HOST');
autoexec;
result=test_true(result,'autoexec direct script dispatch',~isfield(sys,'parallel'),...
                 'direct autoexec script dispatch must leave unknown hosts unchanged');

% Execute the GPU guard pair by script name
spin_system=local_quiet_spin_system();
spin_system.sys.enable={'gpu','mex'};
start_disallow_gpu;
gpu_removed=user_wanted_gpu&&(~ismember('gpu',spin_system.sys.enable));
end_disallow_gpu;
result=test_true(result,'gpu guard direct script dispatch',...
                 gpu_removed&&ismember('gpu',spin_system.sys.enable),...
                 'direct GPU guard script dispatch must remove, and then restore, the gpu flag');

% Execute the parallel profiler pair by script name
spin_system=local_quiet_spin_system();
spin_system.sys.enable={};
spin_system.sys.scratch=tempdir;
parallel_profiler_start;
pause(0.01);
parallel_profiler_report;
result=test_true(result,'parallel profiler direct script dispatch',...
                 exist('nbytes','var')&&exist('walltime','var')&&all(isfinite(nbytes)),...
                 'direct profiler script dispatch must report finite byte counters and wall time');

% Execute the serial Redfield include by script name
[spin_system,Q,L0,R,expected]=local_redfield_fixture(); %#ok<ASGLU>
redfield_integral_serial;
result=test_close(result,'serial Redfield direct script dispatch',R,expected,1e-10,1e-10,...
                  'direct serial Redfield include dispatch must reproduce the analytical integral');

% Execute the asynchronous Redfield include by script name
[spin_system,Q,L0,R,expected]=local_redfield_fixture(); %#ok<ASGLU>
redfield_integral_async;
result=test_close(result,'async Redfield direct script dispatch',R,expected,1e-10,1e-10,...
                  'direct asynchronous Redfield include dispatch must reproduce the analytical integral');

end


function result=local_test_existentials(result,existentials_file)

% Keep a direct reference to the expected production file
result=test_true(result,'existentials canonical path',...
                 strcmp(which('existentials'),existentials_file),...
                 'the read-only existential check must resolve to the canonical Spinach integrity file');

% Run the read-only startup checks from the canonical path setup
startup_text=evalc('local_call_existentials();');
result=test_true(result,'existentials read-only smoke',...
                 contains(startup_text,'Running startup checks'),...
                 'existentials() should complete and announce startup checks on a correctly configured path');

end


function result=local_test_exorcise_patrol(result,exorcise_file,patrol_file)

% Check exorcise input validation without scanning or touching the Wiki
result=test_true(result,'exorcise canonical path',strcmp(which('exorcise'),exorcise_file),...
                 'exorcise must resolve to the production integrity file before validation is tested');
err_text=local_error_text(@()exorcise('bad-mode'));
result=test_true(result,'exorcise invalid mode guard',contains(err_text,'online'),...
                 'exorcise must reject modes other than online and offline before repository scanning starts');

% Check patrol can exit safely on an empty target selection
result=test_true(result,'patrol canonical path',strcmp(which('patrol'),patrol_file),...
                 'patrol must resolve to the production integrity file before validation is tested');
rng_state=rng;
rng_cleanup=onCleanup(@()rng(rng_state));
err_text=local_error_text(@()patrol('__no_such_subject_for_regression_7b31e71d__'));
clear('rng_cleanup');
result=test_true(result,'patrol empty selection exits',contains(err_text,'file list is empty'),...
                 'patrol should fail safely before entering its continuous loop when no example file is selected');

% Verify the high-risk full paths remain intentionally guarded by source
exorcise_src=fileread(exorcise_file);
patrol_src=fileread(patrol_file);
result=test_true(result,'exorcise offline source guard',...
                 contains(exorcise_src,'strcmp(mode,''online'')')&&...
                 contains(exorcise_src,'checkcode(file_name)'),...
                 'the offline mode must skip only the Wiki probe while preserving syntax checking in the source');
result=test_true(result,'patrol continuous-loop source guard',...
                 contains(patrol_src,'while hashes_match')&&contains(patrol_src,'eval(mfiles(n).name'),...
                 'the production patrol owner path is a continuous example runner and is therefore not run in full here');

end


function result=local_test_rearm_sniff_fixture(result,rearm_file,sniff_file)

% Create a synthetic Spinach-like tree under a temporary directory
scratch=tempname(tempdir);
integrity_dir=fullfile(scratch,'kernel','integrity');
mkdir(integrity_dir);
mkdir(fullfile(scratch,'interfaces'));
mkdir(fullfile(scratch,'experiments'));
mkdir(fullfile(scratch,'etc'));
copyfile(rearm_file,fullfile(integrity_dir,'rearm.m'));
copyfile(sniff_file,fullfile(integrity_dir,'sniff.m'));
local_write_dummy_function(fullfile(scratch,'kernel','dummy_integrity_subject.m'));
placeholder=0;
save(fullfile(integrity_dir,'smells.mat'),'placeholder');

% Shadow only inside this fixture and restore afterwards
old_path=path;
old_dir=pwd;
cleanup_obj=onCleanup(@()local_restore_remove(old_path,old_dir,scratch));
addpath(integrity_dir,'-begin');
cd(integrity_dir);
clear('rearm','sniff');

% Rearm the temporary sniffer database without touching the real one
rearm_text=evalc('local_call_rearm();');
loaded=load(fullfile(integrity_dir,'smells.mat'),'smells');
result=test_true(result,'rearm temporary database',...
                 contains(rearm_text,'sniffer rearmed')&&isfield(loaded,'smells')&&...
                 (~isempty(loaded.smells)),...
                 'rearm must overwrite only the fixture smells.mat and populate hash records');

% Sniff the matching temporary tree and check the all-clear path
sniff_text=evalc('sniff(''none'');');
result=test_true(result,'sniff temporary all-clear',contains(sniff_text,'everything smells fine'),...
                 'sniff must accept the freshly rearmed temporary tree');

% Check sniff input validation in the fixture copy
err_text=local_error_text(@()sniff('bad-action'));
result=test_true(result,'sniff invalid action guard',contains(err_text,'none')&&contains(err_text,'open'),...
                 'sniff must reject actions other than none and open');
clear('rearm','sniff');

end


function result=local_test_real_sniff(result,integrity_dir)

% Run the production sniffer read-only from the directory containing smells.mat
old_dir=pwd;
cleanup_obj=onCleanup(@()cd(old_dir));
cd(integrity_dir);
sniff_text=evalc('sniff(''none'');');
result=test_true(result,'sniff production read-only smoke',~isempty(sniff_text),...
                 'production sniff(''none'') must complete read-only and print either all-clear or fishy-file diagnostics');

end


function result=local_test_smack_static(result,smack_file)

% Inspect smack instead of executing its intentionally destructive owner path
smack_handle=@smack; %#ok<NASGU>
smack_src=fileread(smack_file);
result=test_true(result,'smack destructive-path inventory',...
                 contains(smack_src,'delete(gcp(''nocreate''))')&&...
                 contains(smack_src,'parcluster(''Processes'')')&&...
                 contains(smack_src,'fclose(''all'')')&&...
                 contains(smack_src,'clear(''all'')')&&...
                 contains(smack_src,'gpuDeviceCount'),...
                 'smack is command-line cleanup code; static coverage verifies the pool, job, file, workspace, and GPU reset calls without executing them');

end


function local_call_existentials()

% Call the production startup checker directly for output capture
existentials();

end


function local_call_rearm()

% Call the temporary sniffer database rebuilder directly for output capture
rearm();

end


function result=local_test_compile_mex(result,compile_file)

% Create a scratch copy so the no-source branch cannot touch shipped MEX files
scratch=tempname(tempdir);
mkdir(scratch);
copyfile(compile_file,fullfile(scratch,'compile_mex.m'));
old_path=path;
old_dir=pwd;
cleanup_obj=onCleanup(@()local_restore_remove(old_path,old_dir,scratch));
addpath(scratch,'-begin');
clear('compile_mex');

% Exercise the safe early error branch in the scratch copy
err_text=local_error_text(@()compile_mex());
result=test_true(result,'compile_mex no-source guard',contains(err_text,'No C++ source files'),...
                 'compile_mex must refuse to run when its own directory contains no C++ source files');
clear('compile_mex');

% Inspect the production source for platform flags and output confinement
compile_src=fileread(compile_file);
result=test_true(result,'compile_mex source confinement',...
                 contains(compile_src,'dir(fullfile(here,''*.cpp''))')&&...
                 contains(compile_src,'mex(mex_args{:},src_file,''-outdir'',here)')&&...
                 contains(compile_src,'''-R2018a'''),...
                 'compile_mex must discover sources in its own directory, write outputs there, and request the interleaved-complex API');

end


function spin_system=local_quiet_spin_system()

% Build the minimal Spinach object needed by include scripts
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.sys.disable={};
spin_system.sys.scratch=tempdir;

end


function [spin_system,Q,L0,R,expected]=local_redfield_fixture()

% Build the minimal Spinach object needed by Redfield includes
spin_system=local_quiet_spin_system();
spin_system.tols.rlx_integration=1e-6;
spin_system.tols.rlx_zero=1e-14;
spin_system.tols.prop_chop=1e-14;
spin_system.tols.small_matrix=10;
spin_system.tols.dense_matrix=0.5;
spin_system.bas.formalism='sphten-liouv';
spin_system.bas.basis=1;
spin_system.chem.parts={1};
spin_system.rlx.tau_c={1/3};

% Create a one-dimensional rank-one tensor fixture
Q=cell(1,1);
Q{1}=repmat({sparse(1,1)},3,3);
Q{1}{2,2}=sparse(1);
L0=sparse(1,1);
R=sparse(1,1);

% Compute the analytical one-dimensional integral reference
rate=-1;
weight=1/3;
upper_limit=-1.5*(1/rate)*log(1/spin_system.tols.rlx_integration);
expected=-weight*(1-exp(-upper_limit));

end


function created_pool=local_start_pool_if_needed()

% Start a one-worker process pool only when the caller has none
pool=gcp('nocreate');
created_pool=isempty(pool);
if created_pool
    parpool('Processes',1);
end

end


function local_delete_pool(created_pool)

% Delete only the pool created by this draft test
if created_pool
    pool=gcp('nocreate');
    if ~isempty(pool)
        delete(pool);
    end
end

end


function local_restore_autoexec(old_computer,old_position,old_style,old_menu,old_toolbar)

% Restore environment and graphics defaults touched by autoexec
setenv('COMPUTERNAME',old_computer);
set(groot,'defaultFigurePosition',old_position);
set(groot,'defaultFigureWindowStyle',old_style);
set(groot,'defaultFigureMenuBar',old_menu);
set(groot,'defaultFigureToolbar',old_toolbar);

end


function local_write_dummy_function(file_name)

% Write a tiny function into the temporary integrity tree
fid=fopen(file_name,'w');
fprintf(fid,'function answer=dummy_integrity_subject()\n');
fprintf(fid,'answer=1;\n');
fprintf(fid,'end\n\n');
fclose(fid);

end


function local_restore_remove(old_path,old_dir,path_name)

% Restore path and directory before removing the scratch tree
path(old_path);
cd(old_dir);
clear('rearm','sniff','compile_mex');
if exist(path_name,'dir')
    rmdir(path_name,'s');
end

end


function err_text=local_error_text(fun_handle)

% Return an error message without allowing expected failures to escape
try
    fun_handle();
    err_text='';
catch err
    err_text=err.message;
end

end

