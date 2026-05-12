% Tests sparse, tensor-product, and simple graph utilities. Syntax:
%
%                    result=test_sparse_tensor_graph_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks small sparse-format transforms, Kronecker-matrix action,
% graph pruning, permutation group metadata, tuple enumeration, connectivity,
% and simple bin packing against explicit references.
%
% ilya.kuprov@weizmann.ac.il

function result=test_sparse_tensor_graph_suite()

% State the utility target of the test
result=new_test_result('kernel/sparse_tensor_graph_suite',...
                       'Sparse, tensor, and graph helper functions',...
                       'small sparse, tensor-product, graph, and combinatorial helpers must preserve deterministic reference behaviour.');

% Check sparse logical matrix conversion to partial CSR indexing
A=sparse(logical([0 1 0;1 0 1;0 0 0]));
[row_ptr,col_idx]=sparse2csr(A);
result=test_close(result,'sparse2csr row pointer',row_ptr,[1;2;4;4],1e-15,1e-15,...
                  'CSR row pointers are one-based and mark row starts plus the final sentinel');
result=test_close(result,'sparse2csr column index',col_idx,[2;1;3],1e-15,1e-15,...
                  'CSR column indices list each row non-zero in row-major order');

% Check Kronecker-matrix multiplication without opening the product
Q={[1 2;0 -1],[2 0;1 3],[0 1;4 -2]};
X=reshape(1:16,8,2);
K=kron(kron(Q{1},Q{2}),Q{3});
result=test_close(result,'kronm matrix action',kronm(Q,X),K*X,1e-14,1e-14,...
                  'kronm must match explicit Kronecker product multiplication');
result=test_close(result,'kronm_new matrix action',kronm_new(Q,X),K*X,1e-14,1e-14,...
                  'kronm_new must match explicit Kronecker product multiplication');

% Check subgraph pruning removes strict subsets while preserving maximal rows
subgraphs=logical([1 1 0;1 0 0;0 1 1;0 1 0]);
subgraphs_ref=logical([1 1 0;0 1 1]);
result=test_true(result,'prune_subgraphs',isequal(prune_subgraphs(subgraphs),subgraphs_ref),...
                 'subgraphs contained entirely inside larger subgraphs are removed');

% Check a small permutation group table
G=perm_group('S3');
result=test_close(result,'perm_group S3 class sizes',G.class_sizes,[1 2 3],1e-15,1e-15,...
                  'S3 has identity, two 3-cycles, and three transpositions');
result=test_close(result,'perm_group S3 characters',G.class_characters,[1 1 1;1 1 -1;2 -1 0],1e-15,1e-15,...
                  'S3 irreducible character table has the standard three rows');
result=test_true(result,'perm_group S3 order',G.order==sum(G.class_sizes),...
                 'the group order equals the sum of its conjugacy class sizes');

% Check tuple enumeration independent of random output order
rng(1,'twister');
tuples=swizzle({[1 2],[3 4 5]});
tuples_ref=[1 3;1 4;1 5;2 3;2 4;2 5];
result=test_close(result,'swizzle tuple set',sortrows(tuples),tuples_ref,1e-15,1e-15,...
                  'swizzle enumerates the Cartesian product of the supplied index rows');

% Check molecular connectivity on a four-point geometry
xyz=[0 0 0;0.5 0 0;2 0 0;0 0.5 0];
conn_ref=double([0 1 0 1;1 0 0 1;0 0 0 0;1 1 0 0]);
result=test_close(result,'conmat four-point graph',double(conmat(xyz,0.75)),conn_ref,1e-15,1e-15,...
                  'conmat connects points whose Euclidean separation is below the cutoff');

% Check simple first-fit bin packing indices
bins=binpack([4 2 1 5 3],5);
pack_ok=(numel(bins)==4)&&isequal(bins{1},1)&&isequal(bins{2},[2;3])&&...
        isequal(bins{3},4)&&isequal(bins{4},5);
result=test_true(result,'binpack deterministic indices',pack_ok,...
                 'binpack fills bins greedily from the remaining list and returns original indices');

end


