% Radiofrequency-field policy for Kehl ENDOR sequences. Syntax:
%
%      [rf_nutations,rf_auto]=kehl_rf_policy(parameters)
%
% Parameters:
%
%   parameters       - Kehl ENDOR context parameter structure.
%
% Outputs:
%
%   rf_nutations     - user-specified nutation angular-frequency array.
%   rf_auto          - true when nutation fields should be inferred.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_rf_policy.m>

function [rf_nutations,rf_auto]=kehl_rf_policy(parameters)

    % Check consistency
    grumble(parameters);

    % Get explicit nutation angular-frequency data when supplied
    rf_nutations=[];
    if isfield(parameters,'rf_nutations')
        rf_nutations=parameters.rf_nutations;
    end

    % Infer fields from pulses unless explicitly disabled
    if isfield(parameters,'rf_field_from_pulses')
        rf_auto=parameters.rf_field_from_pulses;
    else
        rf_auto=isempty(rf_nutations);
    end

end

% Consistency enforcement
function grumble(parameters)
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
end

