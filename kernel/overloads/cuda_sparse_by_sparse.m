% Sparse matrix product using low-level CUDA CSR arithmetic. Syntax:
%
%                     C=cuda_sparse_by_sparse(A,B,tol)
%
% Parameters:
%
%    A    - real or complex sparse double gpuArray
%
%    B    - real or complex sparse double gpuArray
%
%    tol  - finite non-negative real double scalar; entries with
%           abs(C)<tol are dropped, zero skips dropping
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
% is allocated by a native factory, with columns and numerical
% values written directly into its owned CSR buffers without COO staging.
% Positive tolerances count survivors before allocation and recompute their
% values in bounded shared memory; retained values are not rounded.
% The undocumented R2026b CSR layout is validated before use. A layout
% mismatch is an error, not a fallback to a potentially unsafe GPU product.
% Missing or unloadable platform binaries and GPUs below compute capability
% 7.5 retain native multiplication. Rebuild with the running MATLAB release
% after a CUDA major-version change; the Linux binary uses CUDA 13.
%
% ilya.kuprov@weizmann.ac.il

function C=cuda_sparse_by_sparse(A,B,tol)

% Check consistency
grumble(A,B,tol);

% Return the MATLAB-owned sparse GPU object when the gateway is available
if (exist('cuda_sparse_by_sparse_mex','file')==3)&&...
   (str2double(gpuDevice().ComputeCapability)>=7.5)
    try
        C=cuda_sparse_by_sparse_mex(A,B,tol); return
    catch exception
        if ~strcmp(exception.identifier,'MATLAB:mex:ErrInvalidMEXFile')
            rethrow(exception);
        end
    end
end

% Retain native multiplication and sparse-only dropping without the gateway
C=A*B;
if tol>0
    [rows,cols,values]=find(C);
    keep=~(abs(values)<tol);
    C=sparse(rows(keep),cols(keep),values(keep),size(A,1),size(B,2));
end

end

% Validate the public sparse GPU multiplication interface
function grumble(A,B,tol)
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
if ~isa(tol,'double')||isa(tol,'gpuArray')||~isreal(tol)||...
   ~isscalar(tol)||~isfinite(tol)||(tol<0)
    error('tol must be a finite non-negative real double scalar.');
end
end


