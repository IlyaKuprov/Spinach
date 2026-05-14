% Vacancy-family defect spin systems for diamond. Syntax:
%
%          [sys,inter]=diamond_vacancy(parameters)
%
% Magnetic parameters from Kirui et al., Diam. Relat. Mater. 8,
% 1569 (1999), https://doi.org/10.1016/S0925-9635(99)00037-0
% and Iakoubovskii and Stesmans, Phys. Rev. B 66, 045406 (2002),
% https://doi.org/10.1103/PhysRevB.66.045406; ZFS table values
% cross-checked against Ball, PhD thesis, OIST Graduate University
% (2021).
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .centre       - 'r4_w6', 'w6', 'r4', 'w29', 'r5', 'o1',
%                      'r6', 'r10', or 'r11', default is 'r4_w6'
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
% <https://spindynamics.org/wiki/index.php?title=diamond_vacancy.m>

function [sys,inter]=diamond_vacancy(parameters)

% Set default input
if nargin==0
    parameters=struct();
end

% Check consistency
grumble(parameters);

% Set default parameters
if ~isfield(parameters,'centre')
    parameters.centre='r4_w6';
end
if ~isfield(parameters,'orientation')
    parameters.orientation='111';
end

% Select the vacancy centre
centre=lower(parameters.centre);
switch centre
    case {'r4_w6','w6','r4'}
        electron='E3'; giso=2.0022;
        zfs=diamond_tensor([105 197 -303]*1e6,diamond_frame_z([1 1 1]));
    case 'w29'
        electron='E4'; giso=2.0019;
        zfs=diamond_tensor([297 156 -453]*1e6,diamond_frame_z([1 1 1]));
    case 'r5'
        electron='E3'; giso=2.0023;
        zfs=diamond_tensor([283 244 -524]*1e6,diamond_frame_z([1 1 0]));
    case 'o1'
        electron='E3'; giso=2.0023;
        zfs=diamond_tensor([109 95 -205]*1e6,diamond_frame_z([1 1 0]));
    case 'r6'
        electron='E3'; giso=2.0023;
        zfs=diamond_tensor([62 59 -120]*1e6,diamond_frame_z([1 1 0]));
    case 'r10'
        electron='E3'; giso=2.0023;
        zfs=diamond_tensor([36 36 -73]*1e6,diamond_frame_z([1 1 0]));
    case 'r11'
        electron='E3'; giso=2.0023;
        zfs=diamond_tensor([27 27 -53]*1e6,diamond_frame_z([1 1 0]));
    otherwise
        error('unknown vacancy centre.');
end

% Build the electron tensors
gmat=eye(3)*giso;
nuclei={};

% Build the Spinach structures
[sys,inter]=diamond_system(electron,gmat,zfs,nuclei,parameters.orientation);

end

% Consistency enforcement
function grumble(parameters)
if(~isstruct(parameters))
    error('parameters must be a structure.');
end
if isfield(parameters,'centre')&&(~ischar(parameters.centre))
    error('parameters.centre must be a character string.');
end
if isfield(parameters,'orientation')&&(~ischar(parameters.orientation))
    error('parameters.orientation must be a character string.');
end
end

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

% Enforce symmetric traceless form for quadratic couplings
function M=diamond_traceless(M)
M=(M+M')/2;
M=M-eye(3)*trace(M)/3;
M=(M+M')/2;
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

