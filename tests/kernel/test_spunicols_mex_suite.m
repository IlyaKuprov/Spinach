% Tests the sparse unique-column MEX helper. Syntax:
%
%                    result=test_spunicols_mex_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test compares spunicols() against Matlab unique(A.','rows').' on
% empty, zero-column, duplicate-column, NaN, Inf, signed-value, and random
% sparse real double matrices.
%
% ilya.kuprov@weizmann.ac.il

function result=test_spunicols_mex_suite()

% Announce the test target
fprintf('TESTING: Sparse unique-column MEX helper\n');

% State the utility target of the test
result=new_test_result('kernel/spunicols_mex_suite',...
                       'Sparse unique-column MEX helper',...
                       'spunicols must match Matlab unique(A.'',''rows'').'' for sparse real double matrices.');

% Check empty matrices
A=sparse(0,0); A_ref=unique(A.','rows').';
result=test_true(result,'spunicols empty square',isequal(spunicols(A),A_ref),...
                 'empty sparse matrices should return an empty sparse matrix');
A=sparse(0,5); A_ref=unique(A.','rows').';
result=test_true(result,'spunicols empty wide',isequal(spunicols(A),A_ref),...
                 'zero-row sparse matrices should collapse all empty columns to one column');

% Check zero-column and all-zero matrices
A=sparse(4,0); A_ref=unique(A.','rows').';
result=test_true(result,'spunicols zero columns',isequal(spunicols(A),A_ref),...
                 'zero-column sparse matrices should remain zero-column matrices');
A=sparse(4,3); A_ref=unique(A.','rows').';
result=test_true(result,'spunicols all zero columns',isequal(spunicols(A),A_ref),...
                 'all-zero sparse columns should collapse to one sparse column');

% Check duplicate columns and lexicographic signs
A=sparse([0 0 0 0 0 0;1 1 -1 0 -1 1;0 0 3 0 3 0;-2 -2 0 0 0 -2]);
A_ref=unique(A.','rows').';
result=test_true(result,'spunicols duplicate columns',isequal(spunicols(A),A_ref),...
                 'duplicate sparse columns should match Matlab unique-column ordering');

% Check missing-value ordering and NaN duplicate retention
A=sparse([0 0 NaN NaN 0;Inf Inf 0 0 -Inf;-Inf -Inf 1 1 1;0 0 2 2 0]);
A_ref=unique(A.','rows').';
result=test_true(result,'spunicols nonfinite values',isequaln(spunicols(A),A_ref),...
                 'NaN and Inf ordering must match the Matlab double-transpose reference');

% Check output sparsity
A=sprandn(30,40,0.05); A(:,20:30)=A(:,1:11);
A_mex=spunicols(A);
result=test_true(result,'spunicols sparse output',issparse(A_mex),...
                 'spunicols should return a sparse matrix');

% Check random sparse matrices against Matlab
rng_state=rng;
rng_cleanup=onCleanup(@()rng(rng_state));
rng(1729,'twister');
rand_ok=true;
for n=1:120
    n_rows=randi([1 80]);
    n_cols=randi([1 70]);
    density=10^(-2.7+2.2*rand);
    A=sprandn(n_rows,n_cols,density);
    if (n_cols>1)&&(mod(n,5)==0)
        A(:,n_cols)=A(:,1);
    end
    if (n_cols>2)&&(mod(n,7)==0)
        A(:,2)=sparse(n_rows,1);
    end
    if mod(n,11)==0
        A(randi(n_rows),randi(n_cols))=NaN;
    end
    if mod(n,13)==0
        A(randi(n_rows),randi(n_cols))=Inf;
    end
    A_ref=unique(A.','rows').';
    if ~isequaln(spunicols(A),A_ref)
        rand_ok=false;
        break;
    end
end
result=test_true(result,'spunicols random sparse set',rand_ok,...
                 'random sparse real double matrices should match Matlab unique-column output exactly');
clear('rng_cleanup');

end

