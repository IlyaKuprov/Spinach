% Tests fused multiplication on a selected GPU. Syntax:
%
%                       result=test_multiclean_gpu()
%
% Outputs:
%
%     result  - regression result with value, precision, and
%               device checks; requires an available MATLAB GPU
%
% Sparse GPU panel indexing requires MATLAB R2026b. This test deliberately
% does not select or reset a GPU; unavailable hardware reports SKIP.
%
% ilya.kuprov@weizmann.ac.il

function result=test_multiclean_gpu()

% State the target and explicit local policy
result=new_test_result('kernel/multiclean_gpu','GPU fused multiplication',...
                       'GPU panel arithmetic agrees within roundoff on the execution device.');

% Exercise the supplied global storage policy
spin_system.sys.disable={};
spin_system.tols.small_matrix=200;
spin_system.tols.dense_matrix=0.15;


% Skip unavailable hardware in the default registered test suite
if ~canUseGPU
    result.status='SKIP';
    result.messages{end+1}='No usable GPU is available.';
    return
end

% Cover all storage pairs, both precisions, and real and complex data
for precision={'double','single'}
    for complex_case=0:1
        for left_sparse=0:1
            for right_sparse=0:1
                A=cast(reshape(mod(1:35,7)-3,5,7)/8,precision{1});
                B=cast(reshape(mod(1:77,9)-4,7,11)/8,precision{1});
                if complex_case, A=A+1i*A/2; B=B-1i*B; end
                if left_sparse, A=sparse(A); end
                if right_sparse, B=sparse(B); end
                for device_case=1:3
                    left=A; right=B;
                    if device_case~=2, left=gpuArray(left); end
                    if device_case~=1, right=gpuArray(right); end
                    for tol=[0 1/8 10]
                        label=sprintf('%s complex=%d sparse=%d%d gpu=%d tol=%g',...
                                      precision{1},complex_case,left_sparse,...
                                      right_sparse,device_case,tol);
                        reference=reference_clean(spin_system,left*right,tol);
                        result=compare(result,label,multiclean(spin_system,left,right,tol),reference);
                    end
                end
            end
        end
    end
end

% Preserve sparse output independently of the global product density
A=gpuArray(speye(220));
for column_count=[0 1 33 34 220]
    B=gpuArray(sparse([zeros(220,220-column_count) ones(220,column_count)]));
    result=compare(result,['global density ' num2str(column_count)],...
                   multiclean(spin_system,A,B,1/8),reference_clean(spin_system,A*B,1/8));
end

% Exercise scalar scaling, complex sparse output, and empty shapes
cases={{2,gpuArray(sparse(eye(220)))},...
       {gpuArray(sparse(eye(220))),single(2)},...
       {gpuArray(sparse(complex(eye(220)))),gpuArray(sparse(eye(220)))},...
       {gpuArray(zeros(0,3)),gpuArray(zeros(3,4))},...
       {gpuArray(sparse(4,0)),gpuArray(sparse(0,5))}};
for case_index=1:numel(cases)
    A=cases{case_index}{1}; B=cases{case_index}{2};
    result=compare(result,['GPU shape case ' num2str(case_index)],...
                   multiclean(spin_system,A,B,1/8),reference_clean(spin_system,A*B,1/8));
end

% Disabled clean-up preserves the native GPU reduction
A=gpuArray(randn(31,47)); B=gpuArray(randn(47,53));
result=compare(result,'GPU zero tolerance numerical agreement',multiclean(spin_system,A,B,0),A*B);

% Compare unrounded products within floating-point roundoff
rng(17); A=gpuArray(randn(31,47)); B=gpuArray(randn(47,53));
spin_system.sys.disable={'clean-up'};
result=compare(result,'disabled cleanup numerical agreement',multiclean(spin_system,A,B,1),A*B);
spin_system.sys.disable={};

% Keep a mixed-device left operand resident before sizing panels
for precision={'double','single'}
    A=cast(reshape(1:12,3,4)/8,precision{1});
    B=cast(reshape(1:20,4,5)/8,precision{1});
    result=compare(result,'CPU left with GPU right',...
                   multiclean(spin_system,A,gpuArray(B),1/8),...
                   clean_up(spin_system,A*gpuArray(B),1/8));
    result=compare(result,'GPU left with CPU right',...
                   multiclean(spin_system,gpuArray(A),B,1/8),...
                   clean_up(spin_system,gpuArray(A)*B,1/8));
end

% Global density differs from each panel's local density
A=gpuArray(speye(240)); B=sparse(240,40);
B(:,1:4)=1; B=gpuArray(B);
for density=[0.05 0.15 0.9]
    for small_matrix=[0 200]
        spin_system.tols.dense_matrix=density;
        spin_system.tols.small_matrix=small_matrix;
        A=sparse(A); B=sparse(B);
        for dense_input=0:1
            if dense_input, A=full(A); B=full(B);
            else, A=sparse(A); B=sparse(B); end
            label=sprintf('global density=%g size=%d full=%d',density,small_matrix,dense_input);
            result=compare(result,label,multiclean(spin_system,A,B,0.125),...
                           clean_up(spin_system,A*B,0.125));
        end
    end
end
spin_system.tols.dense_matrix=0.15;
spin_system.tols.small_matrix=200;


end

% Compare GPU results without hiding device transfers or precision changes
function result=compare(result,label,observed,reference)

% Check consistency
grumble(result,label,observed,reference);

% Gather only the small test results for exact comparison
result=test_true(result,[label ' values'],...
                 isequaln(full(gather(observed)),full(gather(reference)))||...
                 gather(norm(observed-reference,'fro'))<=...
                 max(size(reference))*eps(classUnderlying(reference))*...
                 max(1,gather(norm(reference,'fro'))),...
                 'GPU panel clean-up must match the complete product on dyadic inputs.');
result=test_true(result,[label ' precision'],...
                 strcmp(classUnderlying(observed),classUnderlying(reference)),...
                 'The multiplication precision must be preserved.');
result=test_true(result,[label ' device'],isa(observed,'gpuArray'),...
                 'GPU multiplication must return a GPU result.');

end

% Consistency enforcement
function grumble(result,label,observed,reference)
if (~isstruct(result))||(~ischar(label))||(~isrow(label))
    error('a test result and a character row label are required.');
end
if (~isa(observed,'gpuArray'))||(~isa(reference,'gpuArray'))
    error('observed and reference must be GPU arrays.');
end
end


% Reference rounding and global storage policy
function C=reference_clean(spin_system,C,tol)

C=clean_up(spin_system,C,tol);

end


