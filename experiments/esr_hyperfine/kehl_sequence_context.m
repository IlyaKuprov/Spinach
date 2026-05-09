% Sequence-dependent EPR context assembly for Kehl ENDOR. Syntax:
%
%      parameters=kehl_sequence_context(spin_system,parameters)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   parameters       - Kehl ENDOR context parameter structure.
%
% Outputs:
%
%   parameters       - parameter structure with EPR selection data.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_sequence_context.m>

function parameters=kehl_sequence_context(spin_system,parameters)

    % Check consistency
    grumble(spin_system,parameters);

    % Build the EPR sweep axis
    parameters.paramsEPR=kehl_prep_epr(parameters);

    % Select EPR orientations using sequence-specific excitation data
    if parameters.freqDomain==false
        parameters.epr=kehl_ori_field(spin_system,parameters);
    else
        parameters.epr=kehl_ori_freq(spin_system,parameters);
    end

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
    if ~isfield(parameters,'excite_width')
        error('parameters.excite_width must be specified.');
    end
end

