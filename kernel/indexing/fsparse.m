% Fast sparse matrix assembly utility. Syntax:
%
%                    A=fsparse(row_idx,col_idx,vals,n_rows,n_cols)
%
% Parameters:
%
%      row_idx  - row indices, real double or int32 array
%
%      col_idx  - column indices, real double or int32 array
%
%      vals     - real or complex double array, or scalar
%
%      n_rows   - number of matrix rows
%
%      n_cols   - number of matrix columns
%
% Outputs:
%
%      A        - sparse double matrix
%
% This file is a Matlab fallback for the compiled MEX function.
%
% ilya.kuprov@weizmann.ac.il
%
% The MEX implementation uses the sparse-assembly algorithm by
% S. Engblom and D. Lukarski, Parallel Comput. 56, 1-17 (2016).

function A=fsparse(row_idx,col_idx,vals,n_rows,n_cols)

% Check consistency
grumble(row_idx,col_idx,vals,n_rows,n_cols);

% Return Matlab reference sparse matrix
A=sparse(row_idx,col_idx,vals,n_rows,n_cols);

end

% Consistency enforcement
function grumble(row_idx,col_idx,vals,n_rows,n_cols)
if (~isnumeric(row_idx))||(~isreal(row_idx))||issparse(row_idx)||...
   (~ismatrix(row_idx))||(~(isa(row_idx,'double')||isa(row_idx,'int32')))
    error('row_idx must be a real double or int32 matrix.');
end
if (~isnumeric(col_idx))||(~isreal(col_idx))||issparse(col_idx)||...
   (~ismatrix(col_idx))||(~(isa(col_idx,'double')||isa(col_idx,'int32')))
    error('col_idx must be a real double or int32 matrix.');
end
if numel(row_idx)~=numel(col_idx)
    error('row_idx and col_idx must have the same number of elements.');
end
if (~isnumeric(vals))||issparse(vals)||(~isa(vals,'double'))||...
   (~ismatrix(vals))||(~((numel(vals)==numel(row_idx))||...
                         (numel(vals)==1)||...
                         ((numel(vals)==0)&&(numel(row_idx)==0))))
    error('vals must be a full double matrix with one value per index pair, or a scalar.');
end
if (~isnumeric(n_rows))||(~isreal(n_rows))||(~isa(n_rows,'double'))||...
   (~isscalar(n_rows))||(~isfinite(n_rows))||(n_rows<0)||...
   (n_rows~=floor(n_rows))
    error('n_rows must be a non-negative real integer.');
end
if (~isnumeric(n_cols))||(~isreal(n_cols))||(~isa(n_cols,'double'))||...
   (~isscalar(n_cols))||(~isfinite(n_cols))||(n_cols<0)||...
   (n_cols~=floor(n_cols))
    error('n_cols must be a non-negative real integer.');
end
end

