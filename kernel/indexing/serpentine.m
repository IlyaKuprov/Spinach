% Serpentine index matrix used in Spinach for single-index
% numbering of bosonic monomials. Syntax:
%
%                   S=serpentine(nlevels)
%
% Parameters:
%
%    nlevels - number of energy levels in the trun-
%              cated bosonic mode
%
% Outputs:
%
%    idx     - a matrix that indexes the powers of
%              creation and annihilation operators
%              in a grid of bosonic monomials in
%              the following way:
%
%                 (0,0)(0,1)(0,2)     (0)(2)(5)
%                 (1,0)(1,1)(1,2) <=> (1)(4)(7)
%                 (2,0)(2,1)(2,2)     (3)(6)(8)
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=serpentine.m>

function S=serpentine(nlevels)

% Check consistency
grumble(nlevels);

% Build the serpentine index
[rows,cols]=ndgrid(1:nlevels);
[~,idx]=sortrows([rows(:)+cols(:), -rows(:)]);
S=zeros(nlevels,nlevels); S(idx)=0:(nlevels^2-1);

end

% Consistency enforcement
function grumble(nlevels)
if (~isnumeric(nlevels))||(~isreal(nlevels))||...
   (~isscalar(nlevels))||(nlevels<1)||(mod(nlevels,1)~=0)
    error('nlevels must be a positive integer.');
end
end

% There are only two kinds of languages: the ones people
% complain about and the ones nobody uses.
%
% Bjarne Stroustrup,
% creator of C++

