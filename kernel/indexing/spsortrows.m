% Sparse matrix row-sorting permutation utility. Syntax:
%
%                    idx=spsortrows(A)
%
% Parameters:
%
%      A   - sparse real double matrix
%
% Outputs:
%
%      idx - row permutation index, matching the second
%            output of Matlab's sortrows(A)
%
% This file is a Matlab fallback for the compiled MEX function.
%
% ilya.kuprov@weizmann.ac.il

function idx=spsortrows(A)

% Check consistency
grumble(A);

% Return Matlab reference permutation
[~,idx]=sortrows(A);

end

% Consistency enforcement
function grumble(A)
if (~isnumeric(A))||(~issparse(A))||(~isreal(A))||...
   (~isa(A,'double'))||(~ismatrix(A))
    error('A must be a sparse real double matrix.');
end
end


