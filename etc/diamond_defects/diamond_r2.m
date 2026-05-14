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
frame=diamond_frame_z([1 0 0]);
gmat=diamond_tensor([2.0019 2.0019 2.0021],frame);
zfs=diamond_zfs(parameters.d_sign*4173e6,0,frame);
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
if isfield(parameters,'d_sign')&&...
   (~isnumeric(parameters.d_sign)||~isreal(parameters.d_sign)||~isscalar(parameters.d_sign))
    error('parameters.d_sign must be a real scalar.');
end
end

% A self-interstitial has its own symmetry axis.

% Shared local helpers


% Make a principal-axis frame from the z axis
function frame=diamond_frame_z(zaxis)
zaxis=zaxis(:)/norm(zaxis,2);
if abs(dot(zaxis,[0;0;1]))<0.9
    xaxis=cross([0;0;1],zaxis);
else
    xaxis=cross([0;1;0],zaxis);
end
xaxis=xaxis/norm(xaxis,2);
yaxis=cross(zaxis,xaxis);
frame=[xaxis yaxis zaxis];
end

% Orthogonalise a right-handed frame
function frame=diamond_frame_orth(frame)
[frame,~]=qr(frame,0);
if det(frame)<0
    frame(:,3)=-frame(:,3);
end
end

% Build a tensor from principal values and axes
function M=diamond_tensor(values,frame)
frame=diamond_frame_orth(frame);
M=frame*diag(values)*frame';
M=(M+M')/2;
end

% Enforce symmetric traceless form for quadratic couplings
function M=diamond_traceless(M)
M=(M+M')/2;
M=M-eye(3)*trace(M)/3;
M=(M+M')/2;
end

% Build a zero-field splitting tensor from principal axes
function M=diamond_zfs(D,E,frame)
frame=diamond_frame_orth(frame);
M=diamond_traceless(frame*zfs2mat(D,E,0,0,0)*frame');
end

% Crystal-to-laboratory rotation matrix
function C=diamond_orient(orientation)
switch orientation
    case '111'
        C=rotmat_align([1 1 1],[0 0 1]);
    case '110'
        C=rotmat_align([1 1 0],[0 0 1]);
    case '100'
        C=rotmat_align([1 0 0],[0 0 1]);
    otherwise
        error('unknown orientation specification.');
end
end

% Build the Spinach structures
function [sys,inter]=diamond_system(electron,gmat,zfs,nuclei,orientation)
C=diamond_orient(orientation);
sys.isotopes={electron};
inter.zeeman.matrix{1}=C*gmat*C';
if ~isempty(zfs)
    inter.coupling.matrix{1,1}=diamond_traceless(C*zfs*C');
else
    inter.coupling.matrix{1,1}=[];
end
for n=1:numel(nuclei)
    sys.isotopes{n+1}=nuclei{n}.iso;
    inter.zeeman.matrix{n+1}=zeros(3);
    inter.coupling.matrix{1,n+1}=C*nuclei{n}.A*C';
    if ~isempty(nuclei{n}.Q)
        inter.coupling.matrix{n+1,n+1}=diamond_traceless(C*nuclei{n}.Q*C');
    else
        inter.coupling.matrix{n+1,n+1}=[];
    end
end
end

