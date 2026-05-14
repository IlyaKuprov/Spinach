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

% Shared local helpers


% Rotation matrix for an axis-angle rotation in degrees
function R=diamond_rot_axis(axis,angle)
axis=axis(:)/norm(axis,2);
K=[0 -axis(3) axis(2);axis(3) 0 -axis(1);-axis(2) axis(1) 0];
R=eye(3)+sind(angle)*K+(1-cosd(angle))*(K*K);
end

% Orthogonalise a right-handed frame
function frame=diamond_frame_orth(frame)
[frame,~]=qr(frame,0);
if det(frame)<0
    frame(:,3)=-frame(:,3);
end
end

% Make a principal-axis frame from three vectors
function frame=diamond_frame_xyz(xaxis,yaxis,zaxis)
frame=[xaxis(:) yaxis(:) zaxis(:)];
frame=diamond_frame_orth(frame);
end

% Principal-axis frame from a specified x axis and z axis
function frame=diamond_frame_xz(xaxis,zaxis)
xaxis=xaxis(:)/norm(xaxis,2);
zaxis=zaxis(:)-xaxis*dot(xaxis,zaxis(:));
zaxis=zaxis/norm(zaxis,2);
yaxis=cross(zaxis,xaxis);
frame=diamond_frame_xyz(xaxis,yaxis,zaxis);
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

