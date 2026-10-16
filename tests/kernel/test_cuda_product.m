% Compares the CUDA sparse-product MEX with independent CPU products.
% Tests three shapes and all four real/complex operand combinations.
% The tolerance bounds roundoff in both CPU and GPU complex dot products.
% GPU coverage requires MATLAB R2026a or later and a loadable gateway.
% Syntax:
%
%                  result=test_cuda_product()
%
% Outputs:
%
%     result  - regression result with per-case messages and failures
%
% talos@spindynamics.org

function result=test_cuda_product() % #NGRUM

% State the independent reference and report missing GPU prerequisites
result=new_test_result('kernel/cuda_product','CUDA sparse matrix product',...
                       'CUDA sparse products must agree with independent CPU products.');
if isMATLABReleaseOlderThan('R2026a')
    result.messages{end+1}='SKIP: CUDA sparse products require MATLAB R2026a or later.';
    return;
end
if (exist('canUseGPU','file')==0)||(~canUseGPU())
    result.messages{end+1}='SKIP: CUDA sparse products require a usable GPU.';
    return;
end
if str2double(gpuDevice().ComputeCapability)<7
    result.messages{end+1}='SKIP: shipped CUDA gateway requires compute capability 7.0 or later.';
    return;
end
if exist('cuda_sparse_by_sparse_mex','file')~=3
    result.messages{end+1}='SKIP: CUDA sparse products require the compiled gateway.';
    return;
end

% Preserve the caller's random stream and select reproducible inputs
rand_state=rng(); restore_rng=onCleanup(@() rng(rand_state));
rng(19);

% Exercise square, rectangular, and wider sparse products
shapes=[32 32 32;47 61 29;96 80 112];
for n=1:size(shapes,1)
    for k=0:3

        % Generate independent real and imaginary components
        host_a=sprandn(shapes(n,1),shapes(n,2),0.1);
        host_b=sprandn(shapes(n,2),shapes(n,3),0.1);
        if bitget(k,1)
            host_a=host_a+1i*sprandn(shapes(n,1),shapes(n,2),0.1);
        end
        if bitget(k,2)
            host_b=host_b+1i*sprandn(shapes(n,2),shapes(n,3),0.1);
        end

        % Call the gateway directly and form the independent CPU reference
        gpu_a=gpuArray(host_a); gpu_b=gpuArray(host_b);
        try
            gpu_c=cuda_sparse_by_sparse_mex(gpu_a,gpu_b);
        catch exception
            if strcmp(exception.identifier,'MATLAB:mex:ErrInvalidMEXFile')
                result.messages{end+1}='SKIP: CUDA sparse-product gateway cannot be loaded.';
                return;
            end
            rethrow(exception);
        end
        reference=gather(gpu_a)*gather(gpu_b);
        observed=gather(gpu_c);

        % Bound roundoff by eight real operations per complex dot-product term
        roundoff=8*shapes(n,2)*eps;
        gamma=2*roundoff/(1-roundoff);
        tolerance=gamma*norm(abs(host_a)*abs(host_b),'fro');
        label=sprintf('shape %d, complexity %d',n,k);
        result=test_close(result,label,observed,reference,tolerance,0,...
                          'independent sparse CPU and CUDA dot products agree within roundoff');
        result=test_true(result,[label ' storage'],...
                         isa(gpu_c,'gpuArray')&&issparse(gpu_c),...
                         'the gateway returns a sparse gpuArray without a CPU fallback');

    end
end

end


