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
hz_per_mt=abs(spin('E'))/(2*pi)*1e-3;

% Build the electron tensors
electron='E3';
frame=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
            1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
            0          2/sqrt(6)  1/sqrt(3)];
gmat=((frame)*diag([2.0027 2.0027 2.0025])*(frame)');
zfs=frame*zfs2mat(80.3*hz_per_mt,0,0,0,0)*frame';
nuclei={};

% Add the germanium isotope if requested
if isfield(parameters,'germanium')
    germanium=parameters.germanium;
else
    germanium='73Ge';
end
if strcmp(germanium,'73Ge')
    nuclei{end+1}=struct('iso','73Ge','A',eye(3)*1.64*hz_per_mt,'Q',[]);
elseif ~strcmp(germanium,'none')
    nuclei{end+1}=struct('iso',germanium,'A',zeros(3),'Q',[]);
end

% Build the Spinach structures
switch parameters.orientation
    case '111'
        C=rotmat_align([1 1 1],[0 0 1]);
    case '110'
        C=rotmat_align([1 1 0],[0 0 1]);
    case '100'
        C=rotmat_align([1 0 0],[0 0 1]);
    otherwise
        error('unknown orientation specification.');
end
sys.isotopes={electron};
inter.zeeman.matrix{1}=C*gmat*C';
if ~isempty(zfs)
    [~,~,zfs]=mat2ias(C*zfs*C');
    inter.coupling.matrix{1,1}=zfs;
else
    inter.coupling.matrix{1,1}=[];
end
for n=1:numel(nuclei)
    sys.isotopes{n+1}=nuclei{n}.iso;
    inter.zeeman.matrix{n+1}=zeros(3);
    inter.coupling.matrix{1,n+1}=C*nuclei{n}.A*C';
    if ~isempty(nuclei{n}.Q)
        [~,~,nqi]=mat2ias(C*nuclei{n}.Q*C');
        inter.coupling.matrix{n+1,n+1}=nqi;
    else
        inter.coupling.matrix{n+1,n+1}=[];
    end
end

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

