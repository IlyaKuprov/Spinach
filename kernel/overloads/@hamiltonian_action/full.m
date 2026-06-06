% Converts a Hamiltonian action object into a full matrix. Syntax:
%
%                         M=full(H)
%
% Parameters:
%
%     H - a Hamiltonian action object
%
% Outputs:
%
%     M - a full matrix with the same action
%
% ilya.kuprov@weizmann.ac.il
% aditya.dev@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=hamiltonian_action/full.m>

function M=full(H)

% Convert through the sparse representation
M=full(sparse(H));

end

