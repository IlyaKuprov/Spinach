% Nitrogen interstitial spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_n_inter(parameters)
%
% This is a convenience wrapper around diamond_defect().
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .centre        - 'war9' or 'war10', default is 'war9'
%
%      .orientation   - '111', '110', or '100' crystal plane
%                        normal aligned with the magnetic field,
%                        default is '111'
%
%      .include_13c  - include reported 13C hyperfine couplings
%                       where available, false by default
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_n_inter.m>

function [sys,inter]=diamond_n_inter(parameters)

% Set default input
if nargin==0
    parameters=struct();
end

% Check consistency
grumble(parameters);

% Set default centre
if ~isfield(parameters,'centre')
    parameters.centre='war9';
end

% Call the common diamond defect database
parameters.defect=lower(parameters.centre);
[sys,inter]=diamond_defect(parameters);

end

% Consistency enforcement
function grumble(parameters)
if(~isstruct(parameters))
    error('parameters must be a structure.');
end
if isfield(parameters,'centre')&&(~ischar(parameters.centre))
    error('parameters.centre must be a character string.');
end
end

% Good catalogues are maps, not territories.


