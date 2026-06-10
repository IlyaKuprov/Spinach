% Sparse matrix product on the GPU via cuSPARSE ALG3. Syntax:
%
%                 C=cuda_sparse_by_sparse(A,B,chunk_fraction)
%
% Parameters:
%
%    A              - real or complex sparse double gpuArray
%
%    B              - real or complex sparse double gpuArray
%
%    chunk_fraction - cuSPARSE ALG3 chunk fraction, in the range (0,1]
%
% Outputs:
%
%    C              - sparse double gpuArray product A*B
%
% The function uses MATLAB find() to expose GPU-resident sparse triplets,
% calls cuda_sparse_by_sparse_mex() to interpret the column-major sparse
% ordering as transposed CSR storage and run cuSPARSE SpGEMM with
% CUSPARSE_SPGEMM_ALG3, and then reconstructs the sparse gpuArray from the
% resulting triplets.
%
% ilya.kuprov@weizmann.ac.il

function C=cuda_sparse_by_sparse(A,B,chunk_fraction)

% Check consistency
grumble(A,B,chunk_fraction);

% Get matrix dimensions
[n_rows,n_inner]=size(A);
[~,n_cols]=size(B);

% Return an empty GPU sparse matrix when multiplication is vacuous
if (nnz(A)==0)||(nnz(B)==0)
    if isreal(A)&&isreal(B)
        C=sparse(gpuArray.zeros(0,1,'int32'),gpuArray.zeros(0,1,'int32'),...
                 gpuArray.zeros(0,1),n_rows,n_cols);
    else
        C=sparse(gpuArray.zeros(0,1,'int32'),gpuArray.zeros(0,1,'int32'),...
                 complex(gpuArray.zeros(0,1),gpuArray.zeros(0,1)),n_rows,n_cols);
    end
    return
end

% Get sparse triplets on the GPU
[row_a,col_a,val_a]=find(A);
[row_b,col_b,val_b]=find(B);

% Convert one-based GPU indices into zero-based COO arrays
row_a=int32(row_a)-1;
col_a=int32(col_a)-1;
row_b=int32(row_b)-1;
col_b=int32(col_b)-1;

% Pack dimensions for the MEX gateway
dims=uint64([n_rows n_inner n_cols]);

% Run cuSPARSE SpGEMM ALG3 on the GPU
[row_c,col_c,val_c]=cuda_sparse_by_sparse_mex(row_a,col_a,val_a,...
    row_b,col_b,val_b,dims,chunk_fraction);

% Reconstruct the sparse GPU matrix
C=sparse(row_c,col_c,val_c,n_rows,n_cols);

end

function grumble(A,B,chunk_fraction)

if nargin~=3
    error('Three input arguments are required.');
end

if ~isa(A,'gpuArray')
    error('A must be a gpuArray.');
end

if ~isa(B,'gpuArray')
    error('B must be a gpuArray.');
end

if ~issparse(A)
    error('A must be sparse.');
end

if ~issparse(B)
    error('B must be sparse.');
end

if ~strcmp(classUnderlying(A),'double')
    error('A must be double precision.');
end

if ~strcmp(classUnderlying(B),'double')
    error('B must be double precision.');
end

if size(A,2)~=size(B,1)
    error('A and B dimensions are inconsistent.');
end

if ~isfloat(chunk_fraction)||~isreal(chunk_fraction)||...
   ~isscalar(chunk_fraction)||~isfinite(chunk_fraction)
    error('chunk_fraction must be a real finite scalar.');
end

if (chunk_fraction<=0)||(chunk_fraction>1)
    error('chunk_fraction must be in the range (0,1].');
end

max_int=double(intmax('int32'));

if max([size(A,1) size(A,2) size(B,2) nnz(A) nnz(B)])>max_int
    error('Matrix dimensions and nonzero counts must fit into int32.');
end

if exist('cuda_sparse_by_sparse_mex','file')~=3
    error('cuda_sparse_by_sparse_mex is not compiled.');
end

end

