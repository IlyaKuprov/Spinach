% Initial state vector for Kehl ENDOR kernels. Syntax:
%
%      rho0=kehl_rho0(constants,paramsENDOR,B,geff,spin_system,...
%                     parameters,HF_zz,HF_zy,HF_zx,NQI_zz)
%
% Parameters:
%
%   constants        - map containing physical constants.
%   paramsENDOR      - map containing ENDOR parameters.
%   B,geff           - magnetic field and effective g value.
%   spin_system      - Zeeman-Liouville Spinach spin system.
%   parameters       - Kehl ENDOR context parameter structure.
%   HF_zz,HF_zy,HF_zx - effective hyperfine couplings.
%   NQI_zz           - effective quadrupole couplings.
%
% Outputs:
%
%   rho0             - starting Liouville-space state vector.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_rho0.m>

function rho0=kehl_rho0(constants,paramsENDOR,B,geff,spin_system,...
        parameters,HF_zz,HF_zy,HF_zx,NQI_zz)
    if nargin<10
        NQI_zz=0;
    end

    % Check consistency
    grumble(constants,paramsENDOR,B,geff,spin_system,parameters,...
        HF_zz,HF_zy,HF_zx,NQI_zz);

    % Return the cached electron Lz state
    ops=kehl_operator_basis(spin_system,parameters);
    rho0=ops.rho_z;

end

% Consistency enforcement
function grumble(constants,paramsENDOR,B,geff,spin_system,...
        parameters,HF_zz,HF_zy,HF_zx,NQI_zz)
    if ~isa(constants,'containers.Map')
        error('constants must be a containers.Map object.');
    end
    if ~isa(paramsENDOR,'containers.Map')
        error('paramsENDOR must be a containers.Map object.');
    end
    if ~isnumeric(B)
        error('B must be numeric.');
    end
    if ~isnumeric(geff)
        error('geff must be numeric.');
    end
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
    if ~isnumeric(HF_zz)
        error('HF_zz must be numeric.');
    end
    if ~isnumeric(HF_zy)
        error('HF_zy must be numeric.');
    end
    if ~isnumeric(HF_zx)
        error('HF_zx must be numeric.');
    end
    if ~isnumeric(NQI_zz)
        error('NQI_zz must be numeric.');
    end
end

