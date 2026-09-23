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
% Note: integer inputs are evaluated in int64 arithmetic, which is exact
%       for every rank below the square root of intmax('int64'), so that
%       the intermediate L^2+L cannot saturate a narrow integer class or
%       lose precision in double.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=lm2lin.m>

function I=lm2lin(L,M)

% Check consistency
grumble(L,M);

% Get the linear index, exactly in int64 for integer ranks
if isinteger(L)
    I=int64(L).^2+int64(L)-int64(M);
else
    I=L.^2+L-M;
end

% Make sure the index fits the class of the ranks
if isinteger(L)&&any(I(:)>intmax(class(L)))
    error('the linear index does not fit the integer class of L.');
end

% Return in the class of the ranks
I=cast(I,'like',L);

end

% Consistency enforcement
function grumble(L,M)
if (~isnumeric(L))||(~isreal(L))||any(mod(L(:),1)~=0)||...
   (~isnumeric(M))||(~isreal(M))||any(mod(M(:),1)~=0)
    error('all elements of the inputs must be real integers.');
end
if any(abs(M(:))>L(:))
    error('unacceptable projection number.');
end
if any(L(:)<0)
    error('unacceptable total angular momentum.');
end
if any(size(L)~=size(M))
    error('array dimensions are inconsistent.');
end
if isinteger(L)&&any(double(L(:))>=floor(sqrt(double(intmax('int64')))))
    error('integer ranks that large overflow int64 arithmetic.');
end
end

% A casual stroll through the lunatic asylum shows that faith 
% does not prove anything.
%
% Friedrich Nietzsche

