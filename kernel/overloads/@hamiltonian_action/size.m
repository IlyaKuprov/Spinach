% The size of the matrix represented by a Hamiltonian action object
% Syntax:
%
%                         answer=size(H,dim)
%
% Parameters:
%
%    H   - a Hamiltonian action object
%
%    dim - optional, dimension whose size is required
%
% Outputs:
%
%    answer - a vector with one or two elements
%
% ilya.kuprov@weizmann.ac.il
% aditya.dev@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=hamiltonian_action/size.m>

function varargout=size(H,dim)

% Check consistency
if nargin==2, grumble(dim); end

% Compose the answer
if (nargin==1)&&(nargout<=1)
    varargout{1}=[H.dim H.dim];
elseif (nargin==1)&&(nargout==2)
    varargout{1}=H.dim;
    varargout{2}=H.dim;
elseif (nargin==2)&&(dim==1)
    varargout{1}=H.dim;
elseif (nargin==2)&&(dim==2)
    varargout{1}=H.dim;
else
    error('invalid call syntax.');
end

end

% Consistency enforcement
function grumble(dim)
if (~isscalar(dim))||(~ismember(dim,[1 2]))
    error('dim must be 1 or 2.');
end
end

