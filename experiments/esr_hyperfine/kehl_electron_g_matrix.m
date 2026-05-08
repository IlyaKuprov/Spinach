% Electron g tensor for a Spinach spin. Syntax:
%
%      M=kehl_electron_g_matrix(spin_system,spin_idx)
%
% Parameters:
%
%   spin_system      - Spinach spin system structure.
%   spin_idx         - electron spin index.
%
% Outputs:
%
%   M                - electron g tensor.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_electron_g_matrix.m>

function M=kehl_electron_g_matrix(spin_system,spin_idx)

    % Check consistency
    grumble(spin_system,spin_idx);

    % Convert the Spinach Zeeman tensor into a g tensor
    M=spin_system.inter.zeeman.matrix{spin_idx}*...
        spin_system.tols.freeg/spin_system.inter.basefrqs(spin_idx);

end

% Consistency enforcement
function grumble(spin_system,spin_idx)
    if (~isstruct(spin_system))||(~isfield(spin_system,'inter'))||...
            (~isfield(spin_system.inter,'zeeman'))||...
            (~isfield(spin_system.inter.zeeman,'matrix'))||...
            (~isfield(spin_system.inter,'basefrqs'))||(~isfield(spin_system,'tols'))
        error('spin_system must be a Spinach spin system structure.');
    end
    if (~isnumeric(spin_idx))||(~isscalar(spin_idx))||...
            (spin_idx<1)||(mod(spin_idx,1)~=0)
        error('spin_idx must be a positive integer.');
    end
end

