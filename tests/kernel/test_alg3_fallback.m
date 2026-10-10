% Tests optional SpGEMM MEX loading without modifying shipped binaries. Syntax:
%
%                       test_alg3_fallback()
%
% Requires a supported GPU. An isolated copy of the production wrapper
% exercises missing and invalid platform binaries, all four operand
% complexities, and propagation of non-loader errors. No binary is rebuilt.
%
% talos@spindynamics.org

function test_alg3_fallback()

% Isolate the wrapper from every installed SpGEMM binary
saved_path=path; saved_dir=pwd; test_dir=tempname; mkdir(test_dir);
shadow_warning=warning('off','MATLAB:dispatcher:nameConflict');
warning_guard=onCleanup(@()warning(shadow_warning));
wrapper=which('cuda_sparse_by_sparse');
guard=onCleanup(@()restore_test(saved_path,saved_dir,test_dir));
copyfile(wrapper,fullfile(test_dir,'cuda_sparse_by_sparse.m'));
restoredefaultpath; cd(test_dir); clear cuda_sparse_by_sparse;

% Compare native products for each operand complexity without a binary
assert(exist('cuda_sparse_by_sparse_mex','file')==0);
for n=0:3
    A=gpuArray(sparse([1 0 2;0 3 0])*(1+1i*bitget(n,1)));
    B=gpuArray(sparse([1 0;2 3;0 4])*(1+1i*bitget(n,2)));
    C=cuda_sparse_by_sparse(A,B);
    assert(isa(C,'gpuArray')&&issparse(C));
    assert(isequal(gather(C),gather(A*B)));
    assert(isreal(C)==(n==0));
end

% Verify argument validation remains active without the optional binary
rejected=false;
try
    cuda_sparse_by_sparse(full(A),B);
catch exception
    rejected=contains(exception.message,'A must be sparse');
end
assert(rejected);

% Exercise MATLAB's actual invalid-MEX loader with a disposable mock binary
mock_file=fullfile(test_dir,['cuda_sparse_by_sparse_mex.' mexext]);
fid=fopen(mock_file,'w'); fwrite(fid,'not a MEX binary'); fclose(fid);
rehash; assert(exist('cuda_sparse_by_sparse_mex','file')==3);
C=cuda_sparse_by_sparse(A,B);
assert(isequal(gather(C),gather(A*B)));

% Route a mock gateway through the same call boundary without compiling code
delete(mock_file);
fid=fopen(fullfile(test_dir,'exist.m'),'w');
fprintf(fid,['function value=exist(name,kind)\n' ...
             'value=builtin(''exist'',name,kind);\n' ...
             'if strcmp(name,''cuda_sparse_by_sparse_mex''), value=3; end\nend\n']);
fclose(fid);
identifiers={'MATLAB:mex:ErrInvalidMEXFile','Spinach:cuda_sparse_by_sparse_mex:layout',...
             'parallel:gpu:array:OOM','spinach:cuda:compute','spinach:cuda:validation'};
for n=1:numel(identifiers)
    fid=fopen(fullfile(test_dir,'cuda_sparse_by_sparse_mex.m'),'w');
    fprintf(fid,['function value=cuda_sparse_by_sparse_mex(a,b)\n' ...
                 'error(''%s'',''Mock gateway failure.'');\nend\n'],identifiers{n});
    fclose(fid); clear cuda_sparse_by_sparse_mex; rehash;
    if n==1
        C=cuda_sparse_by_sparse(A,B);
        assert(isequal(gather(C),gather(A*B)));
    else
        caught='';
        try
            cuda_sparse_by_sparse(A,B);
        catch exception
            caught=exception.identifier;
        end
        assert(strcmp(caught,identifiers{n}));
    end
end

% Restore the caller environment before reporting completion
clear guard warning_guard;
fprintf('ALG3_FALLBACK_TEST_COMPLETE\n');

end

% Restore path and directory even when an assertion fails
function restore_test(saved_path,saved_dir,test_dir)

% Remove only the disposable test directory and clear its loaded functions
cd(saved_dir); path(saved_path);
clear cuda_sparse_by_sparse cuda_sparse_by_sparse_mex exist;
rmdir(test_dir,'s');

end


