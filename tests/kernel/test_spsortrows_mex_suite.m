% Tests the sparse sortrows MEX helper. Syntax:
%
%                    result=test_spsortrows_mex_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test compares spsortrows() against Matlab sortrows() on empty,
% zero-column, duplicate-row, NaN, Inf, signed-value, and random sparse
% real double matrices.
%
% ilya.kuprov@weizmann.ac.il

function result=test_spsortrows_mex_suite()

% Announce the test target
fprintf('TESTING: Sparse sortrows MEX helper\n');

% State the utility target of the test
result=new_test_result('kernel/spsortrows_mex_suite',...
                       'Sparse sortrows MEX helper',...
                       'spsortrows must return the same row permutation as Matlab sortrows.');

% Check empty matrices
A=sparse(0,0);
[~,idx_ref]=sortrows(A);
result=test_true(result,'spsortrows empty square',isequal(spsortrows(A),idx_ref),...
                 'empty sparse matrices should return an empty permutation');
A=sparse(0,5);
[~,idx_ref]=sortrows(A);
result=test_true(result,'spsortrows empty tall',isequal(spsortrows(A),idx_ref),...
                 'zero-row sparse matrices should preserve Matlab index shape');

% Check zero-column matrices
A=sparse(4,0);
[~,idx_ref]=sortrows(A);
result=test_true(result,'spsortrows zero columns',isequal(spsortrows(A),idx_ref),...
                 'zero-column sparse matrices should keep the input row order');

% Check duplicate rows and lexicographic signs
A=sparse([0 1 0 -2;0 1 0 -2;0 -1 3 0;0 -1 2 9;0 0 0 0]);
[~,idx_ref]=sortrows(A);
result=test_true(result,'spsortrows duplicate rows',isequal(spsortrows(A),idx_ref),...
                 'duplicate sparse rows should keep Matlab stable ordering');

% Check missing-value ordering
A=sparse([0 NaN 0;0 Inf 0;0 -Inf 1;0 0 2;0 NaN -1]);
[~,idx_ref]=sortrows(A);
result=test_true(result,'spsortrows nonfinite values',isequaln(spsortrows(A),idx_ref),...
                 'NaN and Inf ordering must match Matlab sortrows');

% Check random sparse matrices against Matlab
rng_state=rng;
rng_cleanup=onCleanup(@()rng(rng_state));
rng(1729,'twister');
rand_ok=true;
for n=1:80
    n_rows=randi([1 60]);
    n_cols=randi([1 40]);
    density=10^(-2.5+2*rand);
    A=sprandn(n_rows,n_cols,density);
    if mod(n,7)==0
        A(:,1)=sparse(n_rows,1);
    end
    if mod(n,11)==0
        A(randi(n_rows),randi(n_cols))=NaN;
    end
    if mod(n,13)==0
        A(randi(n_rows),randi(n_cols))=Inf;
    end
    [~,idx_ref]=sortrows(A);
    if ~isequal(spsortrows(A),idx_ref)
        rand_ok=false;
        break;
    end
end
result=test_true(result,'spsortrows random sparse set',rand_ok,...
                 'random sparse real double matrices should match Matlab sortrows exactly');
clear('rng_cleanup');

end


