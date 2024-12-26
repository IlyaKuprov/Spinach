% The size of the matrix represented by the OPIUM. Syntax: 
%
%                   answer=size(op,dim)
%
% Parameters:
%
%    op  - an opium object
%
%    dim - optional, dimension whose 
%          size is required
%
% Outputs:
%
%    answer - a vector with one or two elements
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=opium/size.m>

function varargout=size(op,dim)

% Check consistency
if nargin==2, grumble(dim); end

% Compose the answer
if (nargin==1)&&(nargout<=1)
    varargout{1}=[op.dim op.dim];
elseif (nargin==1)&&(nargout==2)
    varargout{1}=op.dim;
    varargout{2}=op.dim;
elseif (nargin==2)&&(dim==1)
    varargout{1}=op.dim;
elseif (nargin==2)&&(dim==2)
    varargout{1}=op.dim;
else
    error('invalid call syntax.');
end

end

% Consistency enforcement
function grumble(dim)
if (~isscalar(dim))||(~ismember(dim,[1 2]))
    error('for a polyadic object, dim must be 1 or 2');
end
end

% Treason doth never prosper: what's the reason?
% Why, if it prosper, none dare call it treason.
%
% John Harrington

