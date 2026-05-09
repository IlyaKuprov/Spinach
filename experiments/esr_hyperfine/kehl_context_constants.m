% Spinach-derived constants for the Kehl ENDOR context. Syntax:
%
%      constants=kehl_context_constants(spin_system)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%
% Outputs:
%
%   constants        - map of physical constants in angular-frequency units.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_context_constants.m>

function constants=kehl_context_constants(spin_system)

    % Check consistency
    grumble(spin_system);

    constants=containers.Map;
    constants('HBAR')=spin_system.tols.hbar;
    constants('K_B')=spin_system.tols.kbol/constants('HBAR');
    constants('MU_B')=spin_system.tols.muB;
    constants('GE')=spin_system.tols.freeg*spin_system.tols.muB/constants('HBAR');
end

% Consistency enforcement
function grumble(spin_system)
    if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
            (~isfield(spin_system,'comp'))
        error('spin_system must be a Spinach spin system structure.');
    end
end

