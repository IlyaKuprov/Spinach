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
% Note: the square root is the only floating-point step; the rest
%       of the arithmetic runs in the class of the input, in the
%       order that cannot overflow. Integer classes are admitted
%       up to the last rank they hold completely, 10 for int8 and
%       180 for int16, the rule by which basis.m picks them.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=lin2lm.m>

function [L,M]=lin2lm(I)

% Check consistency
grumble(I);

% Ranks: integer part of the square root, in the class of the input
L=cast(fix(sqrt(double(I))),'like',I);

% Step back where the floating-point root rounded up to the next integer
rounded_up=(L.^2>I); L(rounded_up)=L(rounded_up)-1;

% Projections, in the evaluation order that cannot overflow
M=L.^2-I+L;

% Make sure every projection lies within its rank
if any(abs(M(:))>L(:))
    error('IEEE arithmetic breakdown, please contact the developer.');
end

end

% Consistency enforcement
function grumble(I)
if (~isnumeric(I))||(~isreal(I))||any(mod(I(:),1)~=0)||any(I(:)<0)
    error('all elements of the input array must be non-negative integers.');
end
if isinteger(I)&&(intmin(class(I))==0)
    error('unsigned integer input is not supported, projections are signed.');
end
if isinteger(I)&&any(I(:)>floor(sqrt(double(intmax(class(I)))+1))^2-1)
    error('integer input beyond the last rank that the class holds completely.');
end
end

% Arrogance on the part of the meritorious is even more offensive
% to us than the arrogance of those without merit: for merit itself
% is offensive.
%
% Friedrich Nietzsche

