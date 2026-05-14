% Field-unit conversion constants for diamond EPR data. Syntax:
%
%          [hz_per_mt,hz_per_t]=diamond_hz_per_mt()
%
% This internal helper returns conversion factors for hyperfine and
% zero-field parameters reported in magnetic-field units.
%
% Outputs:
%
%    hz_per_mt  - free-electron frequency equivalent of one mT, Hz
%
%    hz_per_t   - free-electron frequency equivalent of one T, Hz
%
% <https://spindynamics.org/wiki/index.php?title=diamond_hz_per_mt.m>

function [hz_per_mt,hz_per_t]=diamond_hz_per_mt()

% Convert magnetic-field units using the free electron
hz_per_mt=abs(spin('E'))/(2*pi)*1e-3;
hz_per_t=1e3*hz_per_mt;

end

% Unit conversions should have one home.

