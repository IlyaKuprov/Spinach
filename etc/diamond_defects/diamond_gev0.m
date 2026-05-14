% GeV0 spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_gev0(parameters)
%
% Magnetic parameters from Nadolinny et al., Phys. Status Solidi A
% 213, 2623 (2016), https://doi.org/10.1002/pssa.201600211
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .germanium    - '73Ge', 'none', or another germanium isotope,
%                      default is '73Ge'
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field, default is '111'
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_gev0.m>

function [sys,inter]=diamond_gev0(parameters)

% Set default input
if nargin==0
    parameters=struct();
end

% Check consistency
grumble(parameters);

% Set default orientation
if ~isfield(parameters,'orientation')
    parameters.orientation='111';
end

% Set field-unit conversion constants
[hz_per_mt,~]=diamond_hz_per_mt();

% Build the electron tensors
electron='E3';
frame=diamond_frame_z([1 1 1]);
gmat=diamond_tensor([2.0027 2.0027 2.0025],frame);
zfs=diamond_zfs(80.3*hz_per_mt,0,frame);
nuclei={};

% Add the germanium isotope if requested
germanium=diamond_get(parameters,'germanium','73Ge');
if strcmp(germanium,'73Ge')
    nuclei{end+1}=diamond_nuc('73Ge',eye(3)*1.64*hz_per_mt,[]);
elseif ~strcmp(germanium,'none')
    nuclei{end+1}=diamond_nuc(germanium,zeros(3),[]);
end

% Build the Spinach structures
[sys,inter]=diamond_system(electron,gmat,zfs,nuclei,parameters.orientation);

end

% Consistency enforcement
function grumble(parameters)
if(~isstruct(parameters))
    error('parameters must be a structure.');
end
if isfield(parameters,'orientation')&&(~ischar(parameters.orientation))
    error('parameters.orientation must be a character string.');
end
if isfield(parameters,'germanium')&&(~ischar(parameters.germanium))
    error('parameters.germanium must be a character string.');
end
end

% Germanium follows the group-IV split-vacancy pattern.

