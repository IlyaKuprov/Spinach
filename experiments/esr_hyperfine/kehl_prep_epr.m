% EPR sweep-axis parameters for the Kehl ENDOR context. Syntax:
%
%      paramsEPR=kehl_prep_epr(parameters)
%
% Parameters:
%
%   parameters       - Kehl ENDOR context parameter structure.
%
% Outputs:
%
%   paramsEPR        - map containing EPR field- or frequency-axis data.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_prep_epr.m>

function paramsEPR=kehl_prep_epr(parameters)

    % Check consistency
    grumble(parameters);

    paramsEPR=containers.Map;

    if parameters.freqDomain==true

        % Frequency-domain EPR axis in rad/s
        Npts=round(parameters.epr_freq_range/parameters.epr_freq_step)+1;
        freq=linspace(parameters.epr_freq_min,...
            parameters.epr_freq_min+parameters.epr_freq_range,Npts);
        paramsEPR("freqAxis")=freq;
    else

        % Field-domain EPR axis in T
        fieldCenter=parameters.static_field;
        fieldmin=fieldCenter-0.260;
        fieldmax=fieldCenter+0.260;
        Npts=(round((fieldmax-fieldmin)/parameters.field_step)+1);
        field=linspace(fieldmin,fieldmax,Npts);
        paramsEPR("fieldAxis")=field;
    end

    paramsEPR("Npts")=Npts;

end

% Consistency enforcement
function grumble(parameters)
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
end

