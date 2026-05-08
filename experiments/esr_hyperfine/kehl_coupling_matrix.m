% Spinach coupling tensor lookup converted to Hz. Syntax:
%
%      M=kehl_coupling_matrix(spin_system,row,col)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   row,col          - spin indices in the Spinach coupling matrix.
%
% Outputs:
%
%   M                - coupling tensor in Hz.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_coupling_matrix.m>

function M=kehl_coupling_matrix(spin_system,row,col)

    % Check consistency
    grumble(spin_system,row,col);

    % Get coupling tensor in Hz
    M=kehl_matrix_from_cell(spin_system.inter.coupling.matrix,row,col)/(2*pi);

end

% Consistency enforcement
function grumble(spin_system,row,col)
    if (~isstruct(spin_system))||(~isfield(spin_system,'inter'))||...
            (~isfield(spin_system.inter,'coupling'))||...
            (~isfield(spin_system.inter.coupling,'matrix'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if (~isnumeric(row))||(~isscalar(row))||(row<1)||(mod(row,1)~=0)
        error('row must be a positive integer.');
    end
    if (~isnumeric(col))||(~isscalar(col))||(col<1)||(mod(col,1)~=0)
        error('col must be a positive integer.');
    end
end

