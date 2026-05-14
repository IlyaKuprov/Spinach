% R4/W6 neutral divacancy spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_r4_w6(parameters)
%
% Magnetic parameters from Kirui et al., Diam. Relat. Mater. 8,
% 1569 (1999), https://doi.org/10.1016/S0925-9635(99)00037-0
% ZFS table values cross-checked against Ball, PhD thesis,
% OIST Graduate University (2021).
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field, default is '111'
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_r4_w6.m>

function [sys,inter]=diamond_r4_w6(parameters)

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

% Build the electron tensors
electron='E3';
gmat=eye(3)*2.0022;
zfs=diamond_tensor([105 197 -303]*1e6,diamond_frame_z([1 1 1]));
nuclei={};

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
end

% Vacancy chains leave their spin signature in the D tensor.

