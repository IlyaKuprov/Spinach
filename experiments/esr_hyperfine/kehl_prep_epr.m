% EPR field-axis parameters for the Kehl ENDOR context. Syntax:
%
%      paramsEPR=kehl_prep_epr(spin_system,parameters)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   parameters       - Kehl ENDOR context parameter structure.
%
% Outputs:
%
%   paramsEPR        - map containing EPR field-axis data.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_prep_epr.m>

function paramsEPR=kehl_prep_epr(spin_system,parameters)

    % Check consistency
    grumble(spin_system,parameters);

    % Unpack context data
    constants=parameters.constants;
    g_iso=parameters.g_iso;
    obsField=parameters.field_t;

    gBepr=obsField*g_iso;
    geff=g_iso;

    % EPR simulation range
    fieldCenter=gBepr/g_iso;
    fieldmin=fieldCenter-0.260;
    fieldmax=fieldCenter+0.260;

    % EPR x-axis definition

    % no of points in Field dimension
    Npts=(round((fieldmax-fieldmin)/parameters.field_step_t)+1);
    field=linspace(fieldmin,fieldmax,Npts);

    paramsEPR=containers.Map;
    paramsEPR("Npts")=Npts;
    paramsEPR("fieldAxis")=field;

end
function grumble(spin_system,parameters)
    if (~isstruct(spin_system))||(~isfield(spin_system,'bas'))||(~isfield(spin_system,'comp'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
end

