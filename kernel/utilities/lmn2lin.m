% Converts L,M,N indices of Wigner D functions into linear indices. In 
% the linear indexing convention, Wigner D functions are listed in the
% order of increasing L rank. Within each L, the functions are listed
% in the order of decreasing left index M, and, for each M, in the or-
% der of decreasing N index. One base counting is used:
%
%                     (L=0,M=0,N=0) -> I=1
%                     (L=1,M=1,N=1) -> I=2
%                     (L=1,M=1,N=0) -> I=3, et cetera...
% 
% Syntax:
%
%                         I=lmn2lin(L,M,N)
%
% Parameters:
%
%       L   - ranks of Wigner D functions
%
%       M   - row indices of Wigner D functions
%
%       N   - column indices of Wigner D functions
%
% Outputs:
%
%       I   - linear indices of Wigner D functions, with
%             I=1 corresponding to L=0, M=0, N=0.
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=lmn2lin.m>

function I=lmn2lin(L,M,N)

% Check consistency
grumble(L,M,N);

% Get the linear index
I=L.*(4*L.^2+6*(L-M)+5)/3-M-N+1;

end

% Consistency enforcement
function grumble(L,M,N)
if (~isnumeric(L))||(~isreal(L))||any(mod(L(:),1)~=0)||...
   (~isnumeric(M))||(~isreal(M))||any(mod(M(:),1)~=0)||...
   (~isnumeric(N))||(~isreal(N))||any(mod(N(:),1)~=0)
    error('all elements of the inputs must be real integers.');
end
if any(abs(M(:))>L(:))
    error('unacceptable M projection number.'); 
end
if any(abs(N(:))>L(:))
    error('unacceptable N projection number.'); 
end
if any(L(:)<0)
    error('unacceptable Wigner function rank.'); 
end
if any(size(L)~=size(M))||any(size(L)~=size(N))
    error('array dimensions are inconsistent.');
end
end

% IK has compiled, over the years, a list of books that allow one to
% successfully withstand the toxic social atmosphere of academic in-
% stitutions, and to resist most types of emotional manipulation. In
% the approximate order of reading, the books are:
%
%    - Ayn Rand, "Atlas Shrugged"
%    - Ayn Rand, "The Fountainhead"
%    - Friedrich Nietzsche, "Beyond Good and Evil"
%    - David DeAngelo, "Deep Inner Game"
%    - Ragnar Redbeard, "Might is Right"
%
% If you are starting your career in Academia, or think about apply-
% ing for a faculty post, read these books.

