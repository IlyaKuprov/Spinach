% Converts linear indices of Wigner D functions into L,M,N indices. In 
% the linear indexing convention, Wigner D functions are listed in the
% order of increasing L rank. Within each L, the functions are listed
% in the order of decreasing left index M, and, for each M, in the or-
% der of decreasing N index. One base counting is used:
%
%                     I=1 -> (L=0,M=0,N=0)
%                     I=2 -> (L=1,M=1,N=1)
%                     I=3 -> (L=1,M=1,N=0), et cetera...
% 
% Syntax:
%
%                         [L,M,N]=lin2lmn(I)
%
% Parameters:
%
%       I   - linear indices of Wigner D functions, with
%             I=1 corresponding to L=0, M=0, N=0.
%
% Outputs:
%
%       L   - ranks of Wigner D functions
%
%       M   - row indices of Wigner D functions
%
%       N   - column indices of Wigner D functions
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=lin2lmn.m>

function [L,M,N]=lin2lmn(I)

% Check consistency
grumble(I);

% Solve the cubic equation for L
big_root=(27*I+sqrt(729*I.^2-3)).^(1/3);
L=ceil((3^(1/3)+big_root.^2)./(2*(3^(2/3))*big_root)-1);

% Get the left index
rank_page_position=I-(4*L.^3-L)/3-1;
M=L-fix(rank_page_position./(2*L+1));

% Get the right index
N=L+(2*L+1).*(L-M)-rank_page_position;

% Make sure the conversion is correct
if nnz(lmn2lin(L,M,N)~=I)>0
    error('IEEE arithmetic breakdown, please contact the developer.');
end

end

% Consistency enforcement
function grumble(I)
if (~isnumeric(I))||(~isreal(I))||any(mod(I(:),1)~=0)||any(I(:)<1)
    error('all elements of the input array must be positive integers.');
end
end

% Who gave you the authority to decide which colour 
% Spinach logo was going to be?
%
% Kelly-Anne Ferguson to IK, in April 2011

