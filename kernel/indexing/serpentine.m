% Serpentine index matrix used in Spinach for single-index
% numbering of matrix elements. Syntax:
%
%              S=serpentine(nlevels,idx_base)
%
% Parameters:
%
%    nlevels  - dimension of the matrix, a 
%               positive real integer
% 
%    idx_base - indexing base, 0 or 1 
%
% Outputs:
%
%    S        - serpentine matrix, for example (base 1):
%
%                        (1 )(3 )(6 )(10)
%                        (2 )(5 )(9 )(13)
%                        (4 )(8 )(12)(15)
%                        (7 )(11)(14)(16)
%
%               or (with indexing set to base 0):
%
%                        (0 )(2 )(5 )(9 )
%                        (1 )(4 )(8 )(12)
%                        (3 )(7 )(11)(14)
%                        (6 )(10)(13)(15)
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=serpentine.m>

function S=serpentine(nlevels,idx_base)

% Check consistency
grumble(nlevels,idx_base);

% Build the serpentine matrix
[rows,cols]=ndgrid(1:nlevels);
[~,idx]=sortrows([rows(:)+cols(:), -rows(:)]);
S=zeros(nlevels,nlevels); S(idx)=1:nlevels^2;

% Adjust the indexing base
if idx_base==0, S=S-1; end

end

% Consistency enforcement
function grumble(nlevels,idx_base)
if (~isnumeric(nlevels))||(~isreal(nlevels))||...
   (~isscalar(nlevels))||(nlevels<1)||(mod(nlevels,1)~=0)
    error('nlevels must be a positive integer.');
end
if (~isnumeric(idx_base))||(~isreal(idx_base))||...
   (~isscalar(idx_base))||(~ismember(idx_base,[0 1]))
    error('idx_base must be 0 or 1.');
end
end

% There are only two kinds of languages: the ones people
% complain about and the ones nobody uses.
%
% Bjarne Stroustrup,
% creator of C++

