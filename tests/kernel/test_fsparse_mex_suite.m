% Tests the fast sparse matrix assembly MEX helper. Syntax:
%
%                    result=test_fsparse_mex_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test compares fsparse() against Matlab sparse() on empty, duplicate,
% zero-cancelling, scalar-value, int32-index, complex-value, and random
% triplet inputs.
%
% ilya.kuprov@weizmann.ac.il

function result=test_fsparse_mex_suite()

% Announce the test target
fprintf('TESTING: Fast sparse matrix assembly MEX helper\n');

% State the utility target of the test
result=new_test_result('kernel/fsparse_mex_suite',...
                       'Fast sparse matrix assembly MEX helper',...
                       'fsparse must match Matlab sparse for triplet assembly.');

% Check empty matrices
row_idx=[]; col_idx=[]; vals=[];
A_ref=sparse(row_idx,col_idx,vals,0,5);
result=test_true(result,'fsparse empty wide',isequal(fsparse(row_idx,col_idx,vals,0,5),A_ref),...
                 'empty triplets should assemble an empty sparse matrix of requested size');
A_ref=sparse(row_idx,col_idx,vals,4,0);
result=test_true(result,'fsparse empty tall',isequal(fsparse(row_idx,col_idx,vals,4,0),A_ref),...
                 'empty triplets should support zero-column sparse matrices');

% Check duplicate entries and cancellation
row_idx=[1 1 2 3 3 3 5];
col_idx=[2 2 2 1 1 4 4];
vals=[1 -1 3 5 -2 7 0];
A_ref=sparse(row_idx,col_idx,vals,5,4);
result=test_true(result,'fsparse duplicate cancellation',isequal(fsparse(row_idx,col_idx,vals,5,4),A_ref),...
                 'duplicate indices and zero-cancelling values should match Matlab sparse');

% Check scalar value assembly
row_idx=[1 2 3 3 4];
col_idx=[1 1 2 2 3];
vals=2;
A_ref=sparse(row_idx,col_idx,vals,4,3);
result=test_true(result,'fsparse scalar values',isequal(fsparse(row_idx,col_idx,vals,4,3),A_ref),...
                 'scalar values should expand over all index pairs');

% Check int32 index assembly
row_idx=int32([1 2 2 4 4]);
col_idx=int32([3 1 1 2 2]);
vals=[1 2 3 4 -4];
A_ref=sparse(double(row_idx),double(col_idx),vals,4,3);
result=test_true(result,'fsparse int32 indices',isequal(fsparse(row_idx,col_idx,vals,4,3),A_ref),...
                 'int32 indices should match Matlab sparse with double indices');

% Check complex value fallback
row_idx=[1 2 2 4];
col_idx=[1 1 3 2];
vals=[1+2*1i 3-4*1i -2*1i 5];
A_ref=sparse(row_idx,col_idx,vals,4,3);
result=test_true(result,'fsparse complex values',isequal(fsparse(row_idx,col_idx,vals,4,3),A_ref),...
                 'complex values should match Matlab sparse');

% Check random sparse matrices against Matlab
rng_state=rng;
rng_cleanup=onCleanup(@()rng(rng_state));
rng(1729,'twister');
rand_ok=true;
for n=1:120
    n_rows=randi([1 160]);
    n_cols=randi([1 140]);
    n_vals=randi([1 1200]);
    row_idx=randi(n_rows,n_vals,1);
    col_idx=randi(n_cols,n_vals,1);
    vals=randn(n_vals,1);
    if mod(n,7)==0
        vals(1:min(10,n_vals))=0;
    end
    if mod(n,11)==0
        row_idx=int32(row_idx);
        col_idx=int32(col_idx);
        A_ref=sparse(double(row_idx),double(col_idx),vals,n_rows,n_cols);
    else
        A_ref=sparse(row_idx,col_idx,vals,n_rows,n_cols);
    end
    A=fsparse(row_idx,col_idx,vals,n_rows,n_cols);
    diff_norm=norm(A-A_ref,'fro');
    ref_norm=max(1,norm(A_ref,'fro'));
    tol=100*eps*max(1,n_vals)*ref_norm;
    if (nnz(A)~=nnz(A_ref))||(diff_norm>tol)
        rand_ok=false;
        break;
    end
end
result=test_true(result,'fsparse random sparse set',rand_ok,...
                 'random real double triplet assembly should match Matlab sparse within accumulated roundoff');
clear('rng_cleanup');

end

