% N2V- spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_n2vm(parameters)
%
% Magnetic parameters from Green et al., Phys. Rev. B 92,
% 165204 (2015), https://doi.org/10.1103/PhysRevB.92.165204
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .nitrogen     - '14N' or '15N', default is '15N'
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field, default is '111'
%      .include_13c - include reported 13C hyperfine couplings,
%                     false by default
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_n2vm.m>

function [sys,inter]=diamond_n2vm(parameters)

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
if ~isfield(parameters,'include_13c')
    parameters.include_13c=false;
end

% Build the reported frames
c2rot=diamond_rot_axis([0 0 1],180);
gframe=diamond_frame_xyz([1 1 0],[0 0 1],[1 -1 0]);
nframe=diamond_rot_axis([1 -1 0],-3.5)*...
       diamond_frame_xyz([1 1 -2],[1 1 1],[1 -1 0]);
cframe=diamond_rot_axis([-1 -1 0],2.0)*...
       diamond_frame_xz([-1 -1 0],[-1 1 1]);

% Build the electron tensors
electron='E';
gmat=diamond_tensor([2.00345 2.00274 2.00271],gframe);
zfs=[]; nuclei={};

% Add nitrogen tensors
if strcmp(parameters.nitrogen,'15N')
    nuclei{end+1}=diamond_nuc('15N',diamond_tensor([3.47e6 4.51e6 4.09e6],nframe),[]);
    nuclei{end+1}=diamond_nuc('15N',diamond_tensor([3.47e6 4.51e6 4.09e6],c2rot*nframe),[]);
elseif strcmp(parameters.nitrogen,'14N')
    scale=spin('14N')/spin('15N');
    Qframe=diamond_frame_z([1 1 1]);
    nuclei{end+1}=diamond_nuc('14N',scale*diamond_tensor([3.47e6 4.51e6 4.09e6],nframe),...
                               diamond_zfs(-5.0e6,0,Qframe));
    nuclei{end+1}=diamond_nuc('14N',scale*diamond_tensor([3.47e6 4.51e6 4.09e6],c2rot*nframe),...
                               diamond_zfs(-5.0e6,0,c2rot*Qframe));
else
    error('parameters.nitrogen must be ''14N'' or ''15N''.');
end

% Add reported nearest-neighbour carbons
if parameters.include_13c
    nuclei{end+1}=diamond_nuc('13C',diamond_tensor([202.3e6 202.3e6 317.5e6],cframe),[]);
    nuclei{end+1}=diamond_nuc('13C',diamond_tensor([202.3e6 202.3e6 317.5e6],c2rot*cframe),[]);
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
if isfield(parameters,'nitrogen')&&(~ischar(parameters.nitrogen))
    error('parameters.nitrogen must be a character string.');
end
if isfield(parameters,'include_13c')&&(~islogical(parameters.include_13c))
    error('parameters.include_13c must be logical.');
end
end

% Lower symmetry still deserves a precise frame.

