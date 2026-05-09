% EPR field-axis parameters for the Kehl ENDOR context. Syntax:
%
%      paramsEPR=kehl_prep_epr(parameters)
%
% Parameters:
%
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

function paramsEPR=kehl_prep_epr(parameters)

    % Check consistency
    grumble(parameters);

    % Unpack context data
    g_iso=parameters.g_iso;
    obsField=parameters.field_t;

    gBepr=obsField*g_iso;
    geff=g_iso;

    % EPR simulation range
    fieldCenter=gBepr/geff;
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

% Consistency enforcement
function grumble(parameters)
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
end

