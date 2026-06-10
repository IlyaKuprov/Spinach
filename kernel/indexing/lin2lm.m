% Converts linear indexing of spin states into L,M indexing. In
% the linear indexing convention, spin states are listed in the
% order of increasing L rank, and, within ranks, in the order of
% decreasing M projection. Zero base counting is used: 
%             
%                  I=0 -> (L=0,M=0)
%                  I=1 -> (L=1,M=1)
%                  I=2 -> (L=1,M=0), et cetera...
% 
% Syntax: 
%
%                         [L,M]=lin2lm(I)
%
% Parameters:
%
%       I   - linear indices of spin states, with
%             I=0 corresponding to L=0, M=0.
%
% Outputs:
%
%       L   - ranks of the spin states
%
%       M   - projections of the spin states
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=lin2lm.m>

function [L,M]=lin2lm(I)

% Check consistency
grumble(I);

% A floating point root
L=floor(sqrt(double(I)));

% Back into same type
L=cast(L,'like',I);
M=cast(L.^2+L-I,'like',I);

end

% Consistency enforcement
function grumble(I)
if (~isnumeric(I))||(~isreal(I))||(nnz(mod(I,1))>0)||(nnz(I<0)>0)
    error('the input must contain non-negative real integers.');
end
end

% Arrogance on the part of the meritorious is even more offensive
% to us than the arrogance of those without merit: for merit itself
% is offensive.
%
% Friedrich Nietzsche

