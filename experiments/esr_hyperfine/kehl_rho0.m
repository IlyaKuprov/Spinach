% Initial state vector for Kehl ENDOR kernels. Syntax:
%
%      rho0=kehl_rho0(spin_system,parameters)
%
% Parameters:
%
%   spin_system      - Zeeman-Liouville Spinach spin system.
%   parameters       - Kehl ENDOR context parameter structure.
%
% Outputs:
%
%   rho0             - starting Liouville-space state vector.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_rho0.m>

function rho0=kehl_rho0(spin_system,parameters)

    % Check consistency
    grumble(spin_system,parameters);

    % Request the electron Lz state directly from Spinach
    rho0=-state(spin_system,'Lz',parameters.electron_spin_idx);

end

% Consistency enforcement
function grumble(spin_system,parameters)
    if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||...
            (~isfield(spin_system,'comp'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if ~strcmp(spin_system.bas.formalism,'zeeman-liouv')
        error('spin_system must use zeeman-liouv formalism.');
    end
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
end

