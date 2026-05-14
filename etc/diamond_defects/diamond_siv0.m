% SiV0 spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_siv0(parameters)
%
% Magnetic parameters from Edmonds et al., Phys. Rev. B 77,
% 245205 (2008), https://doi.org/10.1103/PhysRevB.77.245205
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .silicon      - '29Si', 'none', or another silicon isotope,
%                      default is '29Si'
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
% <https://spindynamics.org/wiki/index.php?title=diamond_siv0.m>

function [sys,inter]=diamond_siv0(parameters)

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
if ~isfield(parameters,'include_13c')
    parameters.include_13c=false;
end

% Build the electron tensors
electron='E3';
frame=diamond_frame_z([1 1 1]);
gmat=diamond_tensor([2.0035 2.0035 2.0042],frame);
zfs=diamond_zfs(1000e6,0,frame);
nuclei={};

% Add the silicon isotope if requested
silicon=diamond_get(parameters,'silicon','29Si');
if strcmp(silicon,'29Si')
    nuclei{end+1}=diamond_nuc('29Si',diamond_tensor([78.9e6 78.9e6 76.3e6],frame),[]);
elseif ~strcmp(silicon,'none')
    nuclei{end+1}=diamond_nuc(silicon,zeros(3),[]);
end

% Add reported nearest-neighbour carbons
if parameters.include_13c
    Cmat=diamond_tensor([30.2e6 30.2e6 66.2e6],frame);
    nuc_idx=numel(nuclei);
    nuclei(nuc_idx+1:nuc_idx+6)={[]};
    for n=1:6
        nuclei{nuc_idx+n}=diamond_nuc('13C',Cmat,[]);
    end
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
if isfield(parameters,'silicon')&&(~ischar(parameters.silicon))
    error('parameters.silicon must be a character string.');
end
if isfield(parameters,'include_13c')&&(~islogical(parameters.include_13c))
    error('parameters.include_13c must be logical.');
end
end

% Silicon substitution is visible in both D and hyperfine tensors.

% Shared local helpers


% Field getter with a default value
function value=diamond_get(parameters,name,default)
if isfield(parameters,name)
    value=parameters.(name);
else
    value=default;
end
end

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

% Build a nucleus record
function nucleus=diamond_nuc(iso,A,Q)
nucleus=struct('iso',iso,'A',A,'Q',Q);
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

