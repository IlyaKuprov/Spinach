% Sparse matrix product using low-level CUDA CSR arithmetic. Syntax:
%
%                     C=cuda_sparse_by_sparse(A,B)
%
% Parameters:
%
%    A    - real or complex sparse double gpuArray
%
%    B    - real or complex sparse double gpuArray
%
% Outputs:
%
%    C    - sparse double gpuArray product A*B, complex if either
%           input is complex and neither is all-zero, as in native
%           mtimes; inputs must have compatible sizes
%
% The MEX reads MATLAB's CSR buffers without modifying either input.
% Gustavson row-wise symbolic and numeric phases use bounded shared-memory
% accumulators, with column splitting for overflowing sparse hashes. No
% cuSPARSE multiplication is used. A fresh MATLAB-owned sparse GPU pattern
% is allocated once, and numerical values are written directly into it.
% The undocumented R2026b CSR layout is validated before use. A layout
% mismatch is an error, not a fallback to a potentially unsafe GPU product.
% Missing or unloadable platform binaries retain native multiplication.
%
% ilya.kuprov@weizmann.ac.il

function C=cuda_sparse_by_sparse(A,B)

% Check consistency
grumble(A,B);

% Retain native multiplication when no platform MEX is available
if exist('cuda_sparse_by_sparse_mex','file')~=3
    C=A*B; return
end

% Return the MATLAB-owned sparse GPU object built by the CUDA kernel
try
    C=cuda_sparse_by_sparse_mex(A,B);
catch exception
    if strcmp(exception.identifier,'MATLAB:mex:ErrInvalidMEXFile')
        C=A*B; return
    end
    rethrow(exception);
end

end

% Validate the public sparse GPU multiplication interface
function grumble(A,B)

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

end

