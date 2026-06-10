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
%       L   - ranks of the spin states
%
%       M   - projections of the spin states
%
% Outputs:
%
%       I   - linear indices of spin states, with
%             I=0 corresponding to L=0, M=0.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=lm2lin.m>

function I=lm2lin(L,M)

% Check consistency
grumble(L,M);

% Linear index
I=L.^2+L-M;

end

% Consistency enforcement
function grumble(L,M)
if (~isnumeric(L))||(~isreal(L))||(nnz(mod(L,1))>0)||...
   (~isnumeric(M))||(~isreal(M))||(nnz(mod(M,1))>0)
    error('L and M must contain real integers.');
end
if any(size(L)~=size(M),'all')
    error('dimensions of L and M must be the same.');
end
if nnz(abs(M)>L)>0
    error('each M must be betwen -L and L.');
end
if nnz(L<0)>0
    error('elements of L must be non-negative.');
end
end

% A casual stroll through the lunatic asylum shows that faith 
% does not prove anything.
%
% Friedrich Nietzsche

