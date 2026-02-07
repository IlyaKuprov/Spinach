% Serpentine index matrix used in Spinach for single-index
% numbering of single transition operators. Syntax:
%
%                   S=serpentine(nlevels)
%
% Parameters:
%
%    nlevels - number of energy levels 
%
% Outputs:
%
%    S       - serpentine matrix, for example:
%
%                     (1 )(3 )(6 )(10)
%                     (2 )(5 )(9 )(13)
%                     (4 )(8 )(12)(15)
%                     (7 )(11)(14)(16)
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=serpentine.m>

function S=serpentine(nlevels)

% Check consistency
grumble(nlevels);

% Build the serpentine matrix
[rows,cols]=ndgrid(1:nlevels);
[~,idx]=sortrows([rows(:)+cols(:), -rows(:)]);
S=zeros(nlevels,nlevels); S(idx)=1:nlevels^2;

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

