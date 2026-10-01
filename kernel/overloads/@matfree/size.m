% The size of the matrix represented by the matrix-free core. Syntax:
%
%                   answer=size(op,dim)
%
% Parameters:
%
%    op  - an matfree object
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
% <https://spindynamics.org/wiki/index.php?title=matfree/size.m>

function varargout=size(op,dim)

% Check consistency
if nargin==2, grumble(dim); end

% Compose the answer
if (nargin==1)&&(nargout<=1)
    varargout{1}=op.dims;
elseif (nargin==1)&&(nargout==2)
    varargout{1}=op.dims(1);
    varargout{2}=op.dims(2);
elseif (nargin==2)&&(dim==1)
    varargout{1}=op.dims(1);
elseif (nargin==2)&&(dim==2)
    varargout{1}=op.dims(2);
else
    error('invalid call syntax.');
end

end

% Consistency enforcement
function grumble(dim)
if (~isscalar(dim))||(~ismember(dim,[1 2]))
    error('for a matrix-free core, dim must be 1 or 2');
end
end


