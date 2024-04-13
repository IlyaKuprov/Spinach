% Computes a partial compressed row storage (CSR) transformation 
% for a given Matlab sparse matrix. Adapted from the code written
% by David Gleich. Only returns the index arrays and ignores the
% values. Syntax:
%
%                 [row_ptr,col_idx]=sparse2csr(A)
%
% Parameters:
%
%           A - a Matlab sparse matrix to be converted
%               into the CSR format
%
% Outputs:
%
%           row_ptr - row pointer array of the CSR format
%
%           col_idx - column index array of the CSR format
%
% dgleich@purdue.edu
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=sparse2csr.m>

function [row_ptr,col_idx]=sparse2csr(A)

% Check consistency
grumble(A);

% Set problem dimensions
matrix_dim=size(A,1); n_nonzeros=nnz(A);

% Get Cartesian indices
[rows,cols]=find(A);

% Preallocate the answer
col_idx=zeros(n_nonzeros,1);
row_ptr=zeros(matrix_dim+1,1);

% Count row elements
for n=1:n_nonzeros
    row_ptr(rows(n)+1)=row_ptr(rows(n)+1)+1;
end
row_ptr=cumsum(row_ptr);

% Build column index
for n=1:n_nonzeros
    col_idx(row_ptr(rows(n))+1)=cols(n);
    row_ptr(rows(n))=row_ptr(rows(n))+1;
end

% Build row index
for n=matrix_dim:-1:1
    row_ptr(n+1)=row_ptr(n);
end
row_ptr(1)=0; row_ptr=row_ptr+1;

end

% Consistency enforcement
function grumble(A)
if (~islogical(A))||(~ismatrix(A))||(~issparse(A))
    error('A must be a sparse logical matrix.');
end
end

% It had long since come to my attention that people of
% accomplishment rarely sat back and let things happen
% to them. They went out and happened to things.
%
% Leonardo Da Vinci

