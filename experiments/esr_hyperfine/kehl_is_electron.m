%KEHL_IS_ELECTRON Electron-isotope label test for Kehl ENDOR context.
%
%   Spinach architecture migration May 2026 Talos

function tf=kehl_is_electron(label)
    tf=strcmp(label,'E')||~isempty(regexp(label,'^E\d+$','once'));
end

% Consistency enforcement
