% NV0 excited-state spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_nv0_es(parameters)
%
% Magnetic parameters from Felton et al., Phys. Rev. B 77,
% 081201 (2008), https://doi.org/10.1103/PhysRevB.77.081201
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field, default is '111'
%      .nitrogen     - must be '15N', default is '15N'
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_nv0_es.m>

function [sys,inter]=diamond_nv0_es(parameters)

% Set default input
if nargin==0
    parameters=struct();
end

% Check consistency
grumble(parameters);

% Set default parameters
if ~isfield(parameters,'orientation')
    parameters.orientation='111';
end
if ~isfield(parameters,'nitrogen')
    parameters.nitrogen='15N';
end

% Restrict to the verified isotope
if ~strcmp(parameters.nitrogen,'15N')
    error('NV0 excited-state parameters are only available here for 15N.');
end

% Build the electron tensors
frame=diamond_frame_z([1 1 1]);
electron='E4';
gmat=diamond_tensor([2.0035 2.0035 2.0029],frame);
zfs=diamond_zfs(1685e6,0,frame);

% Add the nitrogen hyperfine tensor
nuclei={diamond_nuc('15N',diamond_tensor([-23.8e6 -23.8e6 -35.7e6],frame),[])};

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
if isfield(parameters,'nitrogen')&&(~ischar(parameters.nitrogen))
    error('parameters.nitrogen must be a character string.');
end
end

% A name is useful only when it remains attached to evidence.

