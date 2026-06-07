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
%    answer - matrix dimensions
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
elseif (nargin==1)
    varargout{1}=H.dim;
    varargout{2}=H.dim;
    for n=3:nargout
        varargout{n}=1;
    end
elseif (nargin==2)
    if dim<=2
        varargout{1}=H.dim;
    else
        varargout{1}=1;
    end
else
    error('invalid call syntax.');
end

end

% Consistency enforcement
function grumble(dim)
if (~isnumeric(dim))||(~isreal(dim))||(~isscalar(dim))||...
   (mod(dim,1)~=0)||(dim<1)
    error('dim must be a positive real integer.');
end
end

