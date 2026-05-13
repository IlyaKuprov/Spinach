% Tests small matrix utility functions. Syntax:
%
%                    result=test_matrix_utility_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks low-level matrix helpers that Spinach uses throughout the
% kernel for algebra, sparsity, block assembly, and indexing operations.
%
% ilya.kuprov@weizmann.ac.il

function result=test_matrix_utility_suite()

% Announce the test target
fprintf('TESTING: Matrix utility functions\n');

% State the utility target of the test
result=new_test_result('kernel/matrix_utility_suite',...
                       'Matrix utility functions',...
                       'matrix helpers must preserve their exact algebraic definitions.');

% Define small test matrices
A=[1 2;3 4];
B=[0 1;-1 2];

% Check anticommutator and cheap norm
result=test_close(result,'acomm',acomm(A,B),A*B+B*A,1e-15,1e-15,...
                  'the anticommutator is AB+BA');
result=test_close(result,'cheap_norm CPU',cheap_norm(A),norm(A,1),1e-15,1e-15,...
                  'for CPU matrices cheap_norm returns the matrix one-norm');

% Check identity and trace predicates
result=test_true(result,'iseye true',iseye(speye(3)),...
                 'the sparse unit matrix is recognised as identity');
result=test_true(result,'iseye false',~iseye([1 0;0 2]),...
                 'a diagonal matrix with a non-unit entry is not identity');
result=test_true(result,'istraceless true',istraceless([1 0;0 -1]),...
                 'zero trace matrices are traceless');
result=test_true(result,'istraceless false',~istraceless(eye(2)),...
                 'the unit matrix has non-zero trace');
result=test_true(result,'krondelta equal',krondelta(3,3),...
                 'Kronecker delta is true for equal integer labels');
result=test_true(result,'krondelta unequal',~krondelta(2,3),...
                 'Kronecker delta is false for unequal integer labels');

% Check sparse block-diagonal assembly
S=sp_block_diag([1 2;3 4],[5;6]);
S_ref=sparse([1 2 0;3 4 0;0 0 5;0 0 6]);
result=test_close(result,'sp_block_diag explicit blocks',S,S_ref,1e-15,1e-15,...
                  'block diagonal assembly places inputs on the diagonal without mixing rows');

% Check row and column replication
M=[1 2 3;4 5 6];
result=test_close(result,'repcols',repcols(M,2,3),[1 2 2 2 3;4 5 5 5 6],1e-15,1e-15,...
                  'repcols repeats selected columns in place');
result=test_close(result,'reprows',reprows(M,1,2),[1 2 3;1 2 3;4 5 6],1e-15,1e-15,...
                  'reprows repeats selected rows in place');

% Check sparse cleanup drops sub-tolerance elements
spin_system.sys.output='hush';
spin_system.sys.disable={};
spin_system.sys.enable={};
spin_system.tols.dense_matrix=0.9;
spin_system.tols.small_matrix=10;
C=sparse([1e-12 1;0 2]);
C_obs=clean_up(spin_system,C,1e-10);
result=test_close(result,'clean_up sparse tolerance',C_obs,sparse([0 1;0 2]),1e-15,1e-15,...
                  'clean_up removes elements below the requested non-zero tolerance');

end

