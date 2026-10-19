% Tests fused multiplication on a selected GPU. Syntax:
%
%                       result=test_multiclean_gpu()
%
% Outputs:
%
%     result  - regression result with value, storage, precision, and
%               device checks; requires an available MATLAB GPU
%
% Sparse GPU panel indexing requires MATLAB R2026b. This test deliberately
% does not select, reset, or silently skip an unavailable GPU.
%
% ilya.kuprov@weizmann.ac.il

function result=test_multiclean_gpu()

% State the target and explicit local policy
result=new_test_result('kernel/multiclean_gpu','GPU fused multiplication',...
                       'GPU panel arithmetic preserves values, storage, and device.');

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
                        reference=reference_clean(left*right,tol);
                        result=compare(result,label,multiclean(left,right,tol),reference);
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
                   multiclean(A,B,1/8),reference_clean(A*B,1/8));
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
                   multiclean(A,B,1/8),reference_clean(A*B,1/8));
end

% Disabled clean-up preserves the native GPU reduction
A=gpuArray(randn(31,47)); B=gpuArray(randn(47,53));
result=compare(result,'GPU zero tolerance native arithmetic',multiclean(A,B,0),A*B);

end

% Compare GPU results without hiding device transfers or precision changes
function result=compare(result,label,observed,reference)

% Check consistency
grumble(result,label,observed,reference);

% Gather only the small test results for exact comparison
result=test_true(result,[label ' values'],...
                 isequaln(full(gather(observed)),full(gather(reference))),...
                 'GPU panel clean-up must match the complete product on dyadic inputs.');
result=test_true(result,[label ' storage'],issparse(observed)==issparse(reference),...
                 'Global representation must not depend on panel dimensions.');
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


% Reference rounding without an application-specific representation policy
function C=reference_clean(C,tol)

% Preserve the storage selected by the native product
spin_system.sys.disable={};
spin_system.tols.small_matrix=0;
if issparse(C)
    spin_system.tols.dense_matrix=Inf;
else
    spin_system.tols.dense_matrix=0;
end
C=clean_up(spin_system,C,tol);

end


