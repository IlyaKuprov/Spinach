% Tests the custom CUDA sparse product against independent CPU products.
% Requires a supported GPU and the compiled gateway. Covers real/complex
% combinations, rectangular shapes, cancellation, stored zeros, wide sparse
% hashes, overflowing rows, MATLAB consumers, and input immutability.
%
% talos@spindynamics.org

function test_cuda_sparse()

% Select the actual production gateway rather than a fallback
assert(exist('cuda_sparse_by_sparse_mex','file')==3);
rng(19);
shapes=[1 1 1 1;5 7 3 0.4;50 40 60 0.1;300 300 300 0.02;...
        0 5 3 0.5;4 0 6 0.5;5 6 4 0;100 1 80 0.5];
for s=1:size(shapes,1)
    for k=0:3
        A=spfun(@(x) round(8*x),sprand(shapes(s,1),shapes(s,2),shapes(s,4)));
        B=spfun(@(x) round(8*x),sprand(shapes(s,2),shapes(s,3),shapes(s,4)));
        if bitget(k,1), A=A+1i*A; end
        if bitget(k,2), B=B-2i*B; end
        Ag=gpuArray(A); Bg=gpuArray(B);
        C=cuda_sparse_by_sparse(Ag,Bg);
        assert(isequal(full(gather(C)),full(A*B)));
        assert(issparse(C)&&isa(C,'gpuArray'));
        assert(isequal(gather(Ag),A)&&isequal(gather(Bg),B));
    end
end

% Check numerical cancellation and explicit GPU storage zeros
A=sparse([1 1 2 2],[1 2 1 2],[1 -1 2 -2],2,2);
B=sparse(ones(2));
C=cuda_sparse_by_sparse(gpuArray(A),gpuArray(B));
assert(nnz(C)==0&&isequal(full(gather(C)),zeros(2)));
Ag=real(gpuArray(complex(A,A)).*gpuArray(sparse([1 0;1 1])));
C=cuda_sparse_by_sparse(Ag,gpuArray(B));
assert(isequal(full(gather(C)),full(gather(Ag))*full(B)));

% Exercise a million-column hash and an overflowing dense row
A=sparse([1 1 2],[1 2 3],[2 3 4],4,3);
B=sparse([1 2 3],[1 500000 1000000],[5 6 7],3,1000000);
C=cuda_sparse_by_sparse(gpuArray(A),gpuArray(B));
assert(isequal(gather(C),A*B));
A=sparse(randi(3,10,30)); B=sparse(randi(3,30,20000));
C=cuda_sparse_by_sparse(gpuArray(A),gpuArray(B));
assert(isequal(gather(C),A*B));

% Check floating-point accuracy and ordinary MATLAB consumers
for k=0:3
    A=sprandn(80,100,0.05)*(1+1i*bitget(k,1));
    B=sprandn(100,80,0.05)*(1+1i*bitget(k,2));
    C=cuda_sparse_by_sparse(gpuArray(A),gpuArray(B));
    R=A*B; tol=1e-13*norm(R,'fro');
    assert(norm(gather(C)-R,'fro')<tol);
    assert(isreal(C)==(k==0));
    assert(norm(gather(C')-R','fro')<tol);
    assert(norm(gather(abs(C))-abs(R),'fro')<tol);
    assert(norm(gather(C*gpuArray.ones(80,1))-R*ones(80,1))<80*tol);
    [row_idx,col_idx,values]=find(C);
    assert(norm(gather(sparse(row_idx,col_idx,values,80,80))-R,'fro')<tol);
    Q=cuda_sparse_by_sparse(C,C);
    assert(norm(gather(Q)-R^2,'fro')<2*norm(R,'fro')*tol);
end

% Check rounding and storage policy using the production cleanup utility
spin_system.sys.disable={};
spin_system.tols.small_matrix=0;
spin_system.tols.dense_matrix=1;
Q=clean_up(spin_system,C,1e-8);
assert(issparse(Q)&&isa(Q,'gpuArray'));
assert(norm(gather(Q)-clean_up(spin_system,R,1e-8),'fro')<tol);

% Exercise stored zeros and allocation capacity after GPU cleanup
A=sprandn(300,300,0.02);
Ag=clean_up(spin_system,gpuArray(A),1);
R=clean_up(spin_system,A,1);
Bg=gpuArray(speye(300));
Q=cuda_sparse_by_sparse(Ag,Bg);
assert(isequal(gather(Q),R));
Q=cuda_sparse_by_sparse(Bg,Ag);
assert(isequal(gather(Q),R));
assert(isequal(gather(Ag),R));

% Reuse an all-zero custom result after native structural compaction
A=sparse([1 1 2 2],[1 2 1 2],[1 -1 2 -2],2,2);
Q=cuda_sparse_by_sparse(gpuArray(A),gpuArray(sparse(ones(2))));
Q=clean_up(spin_system,Q,1e-8);
Q=cuda_sparse_by_sparse(Q,gpuArray(speye(2)));
assert(nnz(Q)==0);
fprintf('CUSTOM_SPARSE_TEST_COMPLETE\n');

end


