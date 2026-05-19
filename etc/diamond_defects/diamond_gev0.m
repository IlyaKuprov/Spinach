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
%      .germanium    - '73Ge', 'none', or another germanium isotope
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_gev0.m>

function [sys,inter]=diamond_gev0(parameters)

% Check input count
if nargin~=1
    error('exactly one input argument is required.');
end

% Check consistency
grumble(parameters);

% Set field-unit conversion constants
hz_per_mt=abs(spin('E'))/(2*pi)*1e-3;

% Build the electron tensors
electron='E3';
frame=[-1/sqrt(2) -1/sqrt(6)  1/sqrt(3);...
        1/sqrt(2) -1/sqrt(6)  1/sqrt(3);...
        0          2/sqrt(6)  1/sqrt(3)];
gmat=((frame)*diag([2.0027 2.0027 2.0025])*(frame)');
zfs=frame*zfs2mat(80.3*hz_per_mt,0,0,0,0)*frame';
nuclei={};

% Add the germanium isotope if requested
germanium=parameters.germanium;
if strcmp(germanium,'73Ge')
    nuclei{end+1}=struct('iso','73Ge','A',eye(3)*1.64*hz_per_mt);
elseif ~strcmp(germanium,'none')
    nuclei{end+1}=struct('iso',germanium);
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
for n=1:numel(nuclei)
    sys.isotopes{n+1}=nuclei{n}.iso;
end
inter.zeeman.matrix=cell(1,numel(sys.isotopes));
inter.zeeman.matrix{1}=C*gmat*C';
inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes));
[~,~,zfs]=mat2ias(C*zfs*C');
inter.coupling.matrix{1,1}=zfs;
for n=1:numel(nuclei)
    if isfield(nuclei{n},'A')&&norm(nuclei{n}.A,2)>0
        inter.coupling.matrix{1,n+1}=C*nuclei{n}.A*C';
    end
end

end

% Consistency enforcement
function grumble(parameters)
if(~isstruct(parameters))
    error('parameters must be a structure.');
end
if ~isfield(parameters,'orientation')
    error('parameters.orientation field is required.');
end
if(~ischar(parameters.orientation))
    error('parameters.orientation must be a character string.');
end
if ~ismember(parameters.orientation,{'111','110','100'})
    error('parameters.orientation must be ''111'', ''110'', or ''100''.');
end
if ~isfield(parameters,'germanium')
    error('parameters.germanium field is required.');
end
if(~ischar(parameters.germanium))
    error('parameters.germanium must be a character string.');
end
end

% The bravest sight in the world is to see a great man 
% struggling against adversity.
% 
% Lucius Annaeus Seneca

