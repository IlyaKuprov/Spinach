% Tests storage-only cleanup after fused matrix multiplication. Syntax:
%
%                    result=test_cleanup_storage()
%
% Outputs:
%
%     result - regression result with rounding and representation checks
%
% ilya.kuprov@weizmann.ac.il

function result=test_cleanup_storage()

% State the test target
result=new_test_result('kernel/cleanup_storage','Exactly-once product cleanup',...
                       'Storage conversion must not round a fused product twice.');

% Distinguish one rounding operation from two at large dynamic range
s.sys.disable={}; s.tols.small_matrix=0; s.tols.dense_matrix=0;
x=-302204600072882.88; tol=0.1;
once=clean_up(s,x,tol); twice=clean_up(s,once,tol);
result=test_true(result,'discriminating rounding case',~isequal(once,twice),...
                 'This input must expose the non-idempotence of floating-point rounding.');
result=test_true(result,'round once',isequal(clean_up(s,once,tol,true),once),...
                 'Storage-only cleanup must retain already-rounded values.');

% Compare complete-product cleanup across all storage combinations
A=diag([1.1 2.2 3.3]); B=diag([0.3 0.7 0.9]);
for sparse_a=[false true]
    for sparse_b=[false true]
        if sparse_a, left=sparse(A); else, left=A; end
        if sparse_b, right=sparse(B); else, right=B; end
        for density=[0 0.5 1]
            s.tols.dense_matrix=density;
            expected=clean_up(s,left*right,tol);
            actual=clean_up(s,multiclean(left,right,tol),tol,true);
            result=test_true(result,'values and global storage',...
                             isequal(actual,expected)&&issparse(actual)==issparse(expected),...
                             'Fused cleanup must retain the whole-matrix storage policy.');
        end
    end
end

% Preserve explicit bypass and recursive storage-only processing
s.sys.disable={'clean-up'};
result=test_true(result,'disabled cleanup',isequal(clean_up(s,x,tol,true),x),...
                 'The disable flag must bypass storage-only cleanup as well.');
s.sys.disable={}; s.tols.dense_matrix=0;
result=test_true(result,'recursive storage-only cleanup',...
                 isequal(clean_up(s,{once,once},tol,true),{once,once}),...
                 'Cell recursion must propagate the already-cleaned flag.');

end

