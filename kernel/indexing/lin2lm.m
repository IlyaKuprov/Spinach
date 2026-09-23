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
%             I=0 corresponding to L=0, M=0; double,
%             single, or a signed integer class
%
% Outputs:
%
%       L   - ranks of the spin states, same class
%             and sparsity as the input
%
%       M   - projections of the spin states, same
%             class and sparsity as the input
%
% Note: the arithmetic runs in double precision so that integer inputs
%       cannot saturate; unsigned integer classes are refused because
%       projections are negative for half of the states.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=lin2lm.m>

function [L,M]=lin2lm(I)

% Check consistency
grumble(I);

% Get the ranks and projections in double precision
L=fix(sqrt(double(I))); M=L.^2+L-double(I);

% Make sure the conversion is correct
if nnz(lm2lin(L,M)~=double(I))>0
    error('IEEE arithmetic breakdown, please contact the developer.');
end

% Return in the class of the input
L=cast(L,'like',I); M=cast(M,'like',I);

end

% Consistency enforcement
function grumble(I)
if (~isnumeric(I))||(~isreal(I))||any(mod(I(:),1)~=0)||any(I(:)<0)
    error('all elements of the input array must be non-negative integers.');
end
if isinteger(I)&&(intmin(class(I))==0)
    error('unsigned integer input is not supported, projections are signed.');
end
if isinteger(I)&&any(I(:)>flintmax)
    error('integer inputs above flintmax are not exactly representable in double precision.');
end
end

% Arrogance on the part of the meritorious is even more offensive
% to us than the arrogance of those without merit: for merit itself
% is offensive.
%
% Friedrich Nietzsche

