% Nuclear chemical-shift tensor for a Spinach spin. Syntax:
%
%      M=kehl_nuclear_cs_matrix(spin_system,spin_idx)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   spin_idx         - nuclear spin index.
%
% Outputs:
%
%   M                - chemical-shift tensor in ppm.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_nuclear_cs_matrix.m>

function M=kehl_nuclear_cs_matrix(spin_system,spin_idx)

    % Check consistency
    grumble(spin_system,spin_idx);

    % Convert the Spinach Zeeman tensor into chemical shift in ppm
    M=(spin_system.inter.zeeman.matrix{spin_idx}/...
        spin_system.inter.basefrqs(spin_idx)-eye(3,3))*1e12;

end

% Consistency enforcement
function grumble(spin_system,spin_idx)
    if (~isstruct(spin_system))||(~isfield(spin_system,'inter'))||...
            (~isfield(spin_system.inter,'zeeman'))||...
            (~isfield(spin_system.inter.zeeman,'matrix'))||...
            (~isfield(spin_system.inter,'basefrqs'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if (~isnumeric(spin_idx))||(~isscalar(spin_idx))||...
            (spin_idx<1)||(mod(spin_idx,1)~=0)
        error('spin_idx must be a positive integer.');
    end
end

