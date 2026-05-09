% Generic field and sweep parameters for the Kehl ENDOR context. Syntax:
%
%      parameters=kehl_context_fields(parameters)
%
% Parameters:
%
%   parameters       - Kehl ENDOR context parameter structure.
%
% Outputs:
%
%   parameters       - parameter structure with checked SI field and sweep data.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_context_fields.m>

function parameters=kehl_context_fields(parameters)

    % Check consistency
    grumble(parameters);

    % Keep caller-supplied SI units without migration

    % Default shaped-pulse multiplicity flag
    if isfield(parameters,'pulse_file')&&~isfield(parameters,'multipulses')
        parameters.multipulses=false;
    end
end

% Consistency enforcement
function grumble(parameters)
    if ~isstruct(parameters)
        error('parameters must be a structure.');
    end
end

