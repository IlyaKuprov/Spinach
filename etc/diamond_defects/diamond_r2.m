% R2 self-interstitial spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_r2(parameters)
%
% Magnetic parameters from Hunt et al., Phys. Rev. B 61,
% 3863 (2000), https://doi.org/10.1103/PhysRevB.61.3863
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .d_sign       - sign of D, default is +1
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field, default is '111'
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_r2.m>

function [sys,inter]=diamond_r2(parameters)

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
if ~isfield(parameters,'d_sign')
    parameters.d_sign=1;
end

% Build the electron tensors
electron='E3';
frame=[0 0 1;...
            1 0 0;...
            0 1 0];
gmat=((frame)*diag([2.0019 2.0019 2.0021])*(frame)');
zfs=frame*zfs2mat(parameters.d_sign*4173e6,0,0,0,0)*frame';
nuclei={};

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
if isfield(parameters,'d_sign')&&...
   (~isnumeric(parameters.d_sign)||~isreal(parameters.d_sign)||~isscalar(parameters.d_sign))
    error('parameters.d_sign must be a real scalar.');
end
end

