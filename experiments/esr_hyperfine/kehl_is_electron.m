% Electron-isotope label test for the Kehl ENDOR context. Syntax:
%
%      tf=kehl_is_electron(label)
%
% Parameters:
%
%   label            - Spinach isotope label.
%
% Outputs:
%
%   tf               - true for electron isotope labels.
%
% February 2024 A. Kehl (akehl@gwdg.de)
% May 2026 Spinach integration
%
% <https://spindynamics.org/wiki/index.php?title=kehl_is_electron.m>

function tf=kehl_is_electron(label)

    % Check consistency
    grumble(label);

    % Detect Spinach electron isotope labels
    tf=strcmp(label,'E')||~isempty(regexp(label,'^E\d+$','once'));

end

% Consistency enforcement
function grumble(label)
    if (~ischar(label))&&(~isstring(label))
        error('label must be a character string or a string scalar.');
    end
end

