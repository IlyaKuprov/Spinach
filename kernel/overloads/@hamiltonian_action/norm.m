% Computes the norm of the matrix represented by a Hamiltonian action
% object. Syntax:
%
%                         n=norm(H,norm_type)
%
% Parameters:
%
%     H         - a Hamiltonian action object
%
%     norm_type - matrix norm type: 1, 2, inf, or 'fro'
%
% Outputs:
%
%     n         - matrix norm
%
% ilya.kuprov@weizmann.ac.il
% aditya.dev@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=hamiltonian_action/norm.m>

function n=norm(H,norm_type) %#NORMOK

% Set the default norm type
if ~exist('norm_type','var'), norm_type=2; end

% Check consistency
grumble(H,norm_type);

% Compute the matrix norm
n=norm(sparse(H),norm_type);

end

% Consistency enforcement
function grumble(H,norm_type)
if ~isa(H,'hamiltonian_action')
    error('H must be a Hamiltonian action object.');
end
if (~(isnumeric(norm_type)&&isreal(norm_type)&&isscalar(norm_type)&&...
      ismember(norm_type,[1 2 inf])))&&...
   (~(ischar(norm_type)&&strcmp(norm_type,'fro')))
    error('norm_type must be 1, 2, inf, or ''fro''.');
end
end

