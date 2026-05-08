%KEHL_RF_POLICY Radiofrequency-field policy for Kehl ENDOR sequences.
%
%   [RF_NUTATIONS,RF_AUTO]=KEHL_RF_POLICY(PARAMETERS) returns explicit
%   nutation frequencies from PARAMETERS.RF_NUTATION_FREQS when supplied,
%   and reports whether the sequence should derive missing nutation fields
%   from pulse durations.
%
%   Inputs:
%
%      PARAMETERS - Kehl ENDOR context parameters.
%
%   Outputs:
%
%      RF_NUTATIONS - user-specified nutation-frequency array.
%      RF_AUTO      - true if nutation fields should be inferred.
%
%   February 2024 A. Kehl (akehl@gwdg.de)
%   Spinach architecture migration May 2026 Talos

function [rf_nutations,rf_auto]=kehl_rf_policy(parameters)

% Check consistency
grumble(parameters);

% Get explicit nutation-frequency data when supplied
rf_nutations=[];
if isfield(parameters,'rf_nutation_freqs')
    rf_nutations=parameters.rf_nutation_freqs;
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
