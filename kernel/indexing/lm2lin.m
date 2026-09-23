% Converts L,M indexing of spin states into linear indexing. In
% the linear indexing convention, spin states are listed in the
% order of increasing L rank, and, within ranks, in the order of
% decreasing M projection. Zero base counting is used: 
%
%                  (L=0,M=0) -> I=0
%                  (L=1,M=1) -> I=1
%                  (L=1,M=0) -> I=2, et cetera...
%
% Syntax: 
%
%                         I=lm2lin(L,M)
%
% Parameters:
%
%       L   - ranks of the spin states; double, single,
%             or a signed integer class
%
%       M   - projections of the spin states, same class
%             as L
%
% Outputs:
%
%       I   - linear indices of spin states, with I=0
%             corresponding to L=0, M=0; same class and
%             sparsity as L
%
% Note: the arithmetic runs in the class of the inputs, in the order
%       that cannot overflow an integer class holding every index of
%       the highest rank present; the grumbler enforces that bound.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=lm2lin.m>

function I=lm2lin(L,M)

% Check consistency
grumble(L,M);

% Linear index, in the evaluation order that cannot overflow
I=L.^2-M+L;

end

% Consistency enforcement
function grumble(L,M)
if (~isnumeric(L))||(~isreal(L))||any(mod(L(:),1)~=0)||...
   (~isnumeric(M))||(~isreal(M))||any(mod(M(:),1)~=0)
    error('all elements of the inputs must be real integers.');
end
if ~strcmp(class(L),class(M))
    error('L and M must have the same class.');
end
if any(size(L)~=size(M))
    error('array dimensions are inconsistent.');
end
if any(abs(M(:))>L(:))
    error('unacceptable projection number.');
end
if any(L(:)<0)
    error('unacceptable total angular momentum.');
end
if isinteger(L)&&any(L(:)>floor(sqrt(double(intmax(class(L)))+1))-1)
    error('the integer class of L cannot hold the indices of its highest rank.');
end
end

% A casual stroll through the lunatic asylum shows that faith 
% does not prove anything.
%
% Friedrich Nietzsche

