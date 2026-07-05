% Sparse matrix unique-column utility. Syntax:
%
%                    A=spunicols(A)
%
% Parameters:
%
%      A   - sparse real double matrix
%
% Outputs:
%
%      A   - sparse real double matrix containing
%            unique columns of the input matrix
%
% This file is a Matlab fallback for the compiled MEX function.
%
% ilya.kuprov@weizmann.ac.il

function A=spunicols(A)

% Check consistency
grumble(A);

% Return Matlab reference unique columns
A=unique(A.','rows').';

end

% Consistency enforcement
function grumble(A)
if (~isnumeric(A))||(~issparse(A))||(~isreal(A))||...
   (~isa(A,'double'))||(~ismatrix(A))
    error('A must be a sparse real double matrix.');
end
end

