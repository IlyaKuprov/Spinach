% Tests fused panel multiplication against the existing clean_up. Syntax:
%
%                         result=test_multiclean()
%
% Outputs:
%
%     result  - regression result with per-case numerical, storage,
%               precision, and shape checks
%
% Uses exactly representable dyadic inputs to isolate clean-up semantics
% from legitimate last-bit differences between blocked BLAS products.
%
% ilya.kuprov@weizmann.ac.il

function result=test_multiclean()

% State the test target
result=new_test_result('kernel/multiclean','Fused multiplication and clean-up',...
                       'Panel products preserve clean_up values and native sparse/full storage.');


% Exercise every sparse/full combination in both floating precisions
A=reshape(mod(1:35,7)-3,5,7)/8;
B=reshape(mod(1:77,9)-4,7,11)/8;
for precision={'double','single'}
    for complex_case=0:1
        for left_sparse=0:1
            for right_sparse=0:1
                left=cast(A,precision{1}); right=cast(B,precision{1});
                if complex_case, left=left+1i*left/2; right=right-1i*right; end
                if left_sparse, left=sparse(left); end
                if right_sparse, right=sparse(right); end
                for tol=[0 1/8 10]
                    label=sprintf('%s complex=%d sparse=%d%d tol=%g',...
                                  precision{1},complex_case,left_sparse,right_sparse,tol);
                    reference=reference_clean(left*right,tol);
                    result=compare(result,label,multiclean(left,right,tol),reference);
                end
            end
        end
    end
end

% Exercise density, dimension, and zero-matrix representation boundaries
for dimension=[199 200 201]
    for occupied=[0 1 30 31 201]
        column_count=min(occupied,dimension);
        A=speye(dimension);
        B=sparse([ones(dimension,column_count) zeros(dimension,201-column_count)]);
        for dense_input=0:1
            if dense_input, A=full(A); B=full(B); end
            label=sprintf('size=%d occupied=%d full=%d',dimension,occupied,dense_input);
            reference=reference_clean(A*B,1/8);
            result=compare(result,label,multiclean(A,B,1/8),reference);
        end
    end
end

% Exercise exactly the 0.15 density boundary on a sufficiently large matrix
A=speye(220); B=sparse([ones(220,33) zeros(220,187)]);
result=compare(result,'sparse exact density',multiclean(A,B,1/8),...
               reference_clean(A*B,1/8));
result=compare(result,'full exact density',multiclean(full(A),full(B),1/8),...
               reference_clean(full(A)*full(B),1/8));

% Clean-up rounds components, including ties, rather than dropping by magnitude
A=[-3/16 -1/16 1/16 3/16 1/32+1i/16];
result=compare(result,'rounding ties and complex components',multiclean(1,A,1/8),...
               reference_clean(A,1/8));

% Repeating clean_up is not generally idempotent in floating point
A=-302204600072882.88;
result=compare(result,'exactly one rounding pass',multiclean(A,1,0.1),...
               reference_clean(A,0.1));

% Cancellation is completed before any rounding is applied
A=[1 1]; B=[3/64 3/64;3/64 -3/64];
result=compare(result,'complete inner products',multiclean(A,B,1/8),...
               reference_clean(A*B,1/8));

% Exercise scalar multiplication, vectors, and zero inner dimensions
cases={{2,sparse(A)},{sparse(A),2},{single(2),single(B)},...
       {zeros(0,3),zeros(3,4)},{zeros(4,3),zeros(3,0)},...
       {zeros(4,0),zeros(0,5)},{sparse(4,0),sparse(0,5)},...
       {zeros(0,0),2},{2,zeros(0,0)},{(1:4)',1:7},...
       {single(eye(3)),double(eye(3))},...
       {sparse([1 0 1]),sparse(reshape(mod(1:900,7),3,300))},...
       {sparse(reshape(mod(1:900,7),300,3)),sparse([1;0;1])}};
for case_index=1:numel(cases)
    A=cases{case_index}{1}; B=cases{case_index}{2};
    result=compare(result,['shape case ' num2str(case_index)],...
                   multiclean(A,B,1/8),reference_clean(A*B,1/8));
end

% A single tolerance can change precision through existing clean_up arithmetic
for sparse_case=0:1
    A=eye(220); if sparse_case, A=sparse(A); end
    for tol=single([0 1/8])
        result=compare(result,'single tolerance precision',multiclean(A,A,tol),...
                       reference_clean(A*A,tol));
    end
end

% Preserve nonfinite arithmetic rather than filtering exceptional entries
A=[Inf NaN;1 -Inf]; B=eye(2);
result=compare(result,'nonfinite entries',multiclean(A,B,1/8),...
               reference_clean(A*B,1/8));

% Large non-dyadic entries must be rounded only once, including sparse panels
for sparse_case=0:1
    A=[-302204600072882.88 302204600072882.88];
    if sparse_case, A=sparse(A); end
    once=reference_clean(A,0.1); twice=reference_clean(once,0.1);
    result=test_true(result,'adversarial rounding is non-idempotent',...
                     ~isequaln(once,twice),'The case must distinguish repeated rounding.');
    result=compare(result,'single rounding of large entries',...
                   multiclean(1,A,0.1),once);
end

% Disabled clean-up returns native BLAS arithmetic, not panel arithmetic
rng(17); A=randn(31,47); B=randn(47,53);
result=compare(result,'zero tolerance native arithmetic',multiclean(A,B,0),A*B);

% Disabled clean-up returns native BLAS arithmetic, not panel arithmetic
rng(17); A=randn(31,47); B=randn(47,53);
result=compare(result,'zero tolerance native arithmetic',multiclean(A,B,0),A*B);

% Reject malformed inputs without accepting ambiguous tolerance policies
bad_cases={{ones(2),ones(3),1},{ones(2),ones(2),-1},...
           {ones(2),ones(2),NaN},{ones(2),ones(2),Inf},...
           {ones(2),ones(2),1i},{ones(2),ones(2),[1 2]},...
           {int32(1),1,1},{ones(2,2,2),1,1},{1,1,int32(1)}};
for case_index=1:numel(bad_cases)
    rejected=false;
    try
        multiclean(bad_cases{case_index}{:});
    catch
        rejected=true;
    end
    result=test_true(result,['invalid case ' num2str(case_index)],rejected,...
                     'Invalid arguments must be rejected.');
end

end

% Compare numerical values and global representation without tolerance masking
function result=compare(result,label,observed,reference)

% Check consistency
grumble(result,label,observed,reference);

% Exact dyadic examples separate rounding errors from BLAS summation order
result=test_true(result,[label ' values'],isequaln(full(observed),full(reference)),...
                 'Panel clean-up must agree with clean_up of the complete product.');
result=test_true(result,[label ' storage'],issparse(observed)==issparse(reference),...
                 'Native product sparse/full representation must be preserved.');
result=test_true(result,[label ' precision'],strcmp(class(observed),class(reference)),...
                 'The multiplication precision must be preserved.');

end

% Consistency enforcement
function grumble(result,label,observed,reference)
if (~isstruct(result))||(~ischar(label))||(~isrow(label))
    error('a test result and a character row label are required.');
end
if (~isnumeric(observed))||(~isnumeric(reference))
    error('observed and reference must be numeric arrays.');
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


