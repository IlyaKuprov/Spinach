% Sparse matrix product on the GPU via cuSPARSE SpGEMM. Syntax:
%
%                     C=cuda_sparse_by_sparse(A,B,alg)
%
% Parameters:
%
%    A    - real or complex sparse double gpuArray
%
%    B    - real or complex sparse double gpuArray
%
%    alg  - cuSPARSE SpGEMM algorithm, 1, 2, or 3 for
%           CUSPARSE_SPGEMM_ALG1, ALG2, or ALG3
%
% Outputs:
%
%    C    - sparse double gpuArray product A*B, complex if either
%           input is complex and neither is all-zero, as in native
%           mtimes; inputs must have compatible sizes
%
% The function passes A and B to cuda_sparse_by_sparse_mex(), which
% reads MATLAB's internal CSR storage of the sparse gpuArrays in place,
% runs cuSPARSE SpGEMM, and returns the product as row-major triplets
% from which the sparse gpuArray is assembled. If the platform MEX is
% absent or MATLAB cannot load it, or if the MEX does not recognise the
% internal storage layout of this MATLAB version, native GPU multiplica-
% tion is used. Other failures are not intercepted.
%
% ilya.kuprov@weizmann.ac.il

function C=cuda_sparse_by_sparse(A,B,alg)

% Check consistency
grumble(A,B,alg);

% Retain native multiplication when no platform MEX is available
if exist('cuda_sparse_by_sparse_mex','file')~=3
    C=A*B; return
end

% Run cuSPARSE SpGEMM on the GPU
try
    [row_c,col_c,val_c]=cuda_sparse_by_sparse_mex(A,B,alg);
catch exception
    if ismember(exception.identifier,{'MATLAB:mex:ErrInvalidMEXFile',...
                                      'Spinach:cuda_sparse_by_sparse_mex:layout'})
        C=A*B; return
    end
    rethrow(exception);
end

% Assemble the sparse GPU matrix
C=sparse(row_c,col_c,val_c,size(A,1),size(B,2));

end

% Validate the public sparse GPU multiplication interface
function grumble(A,B,alg)

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

if ~isa(alg,'double')||isa(alg,'gpuArray')||~isscalar(alg)||~ismember(alg,[1 2 3])
    error('alg must be a CPU double scalar equal to 1, 2, or 3.');
end

end


