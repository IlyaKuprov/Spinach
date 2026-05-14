% Phosphorus-related defect spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_p(parameters)
%
% Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),
% https://doi.org/10.3390/cryst7080237
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .centre       - 'ma1', 'np1', 'np2', 'np3', 'np4', 'np5',
%                      'np6', 'np8', or 'np9', default is 'ma1'
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field, default is '111'
%      .include_13c - include reported 13C hyperfine couplings for MA1,
%                     false by default
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_p.m>

function [sys,inter]=diamond_p(parameters)

% Set default input
if nargin==0
    parameters=struct();
end

% Check consistency
grumble(parameters);

% Set default parameters
if ~isfield(parameters,'centre')
    parameters.centre='ma1';
end
if ~isfield(parameters,'orientation')
    parameters.orientation='111';
end
if ~isfield(parameters,'include_13c')
    parameters.include_13c=false;
end

% Set field-unit conversion constants
[hz_per_mt,~]=diamond_hz_per_mt();

% Select the phosphorus centre
electron='E'; zfs=[];
[gmat,nuclei]=phosphorus_data(lower(parameters.centre),hz_per_mt,parameters.include_13c);

% Build the Spinach structures
[sys,inter]=diamond_system(electron,gmat,zfs,nuclei,parameters.orientation);

end

% Phosphorus centre table data
function [gmat,nuclei]=phosphorus_data(centre,hz_per_mt,include_13c)
nuclei={}; frame=eye(3);
switch centre
    case 'ma1'
        gmat=eye(3)*2.0025;
        nuclei{end+1}=diamond_nuc('31P',...
            diamond_tensor([1.96 1.96 2.32]*hz_per_mt,diamond_frame_z([1 1 1])),[]);
        if include_13c
            nuclei{end+1}=diamond_nuc('13C',...
                diamond_tensor([13.92 13.92 18.13]*hz_per_mt,diamond_frame_z([1 1 1])),[]);
        end
    case 'np1'
        gmat=diamond_tensor([2.00243 2.0028 2.0026],frame);
        nuclei{end+1}=diamond_nuc('31P',diamond_tensor([2.08 2.02 2.18]*hz_per_mt,frame),[]);
        nuclei{end+1}=diamond_nuc('14N',diamond_tensor([4.08 3.10 3.00]*hz_per_mt,frame),[]);
    case 'np2'
        gmat=eye(3)*2.0025;
        nuclei{end+1}=diamond_nuc('31P',...
            diamond_tensor([2.09 2.09 2.34]*hz_per_mt,diamond_frame_z([1 1 1])),[]);
        nuclei{end+1}=diamond_nuc('14N',...
            diamond_tensor([3.09 3.09 6.42]*hz_per_mt,diamond_frame_z([1 1 1])),[]);
    case 'np3'
        gmat=eye(3)*2.0025;
        nuclei{end+1}=diamond_nuc('31P',...
            diamond_tensor([18.23 18.23 17.48]*hz_per_mt,diamond_frame_z([1 1 1])),[]);
        nuclei{end+1}=diamond_nuc('14N',...
            diamond_tensor([0.33 0.33 0.10]*hz_per_mt,diamond_frame_z([1 1 1])),[]);
    case 'np4'
        gmat=diamond_tensor([2.0009 2.0012 2.00047],frame);
        nuclei{end+1}=diamond_nuc('31P',diamond_tensor([5.456 3.838 3.80]*hz_per_mt,frame),[]);
    case 'np5'
        gmat=diamond_tensor([2.0009 2.0009 2.00087],diamond_frame_z([1 1 1]));
        nuclei{end+1}=diamond_nuc('31P',...
            diamond_tensor([1.024 1.024 6.522]*hz_per_mt,diamond_frame_z([1 1 1])),[]);
    case 'np6'
        gmat=diamond_tensor([2.00083 2.00083 2.00085],diamond_frame_z([1 1 1]));
        nuclei{end+1}=diamond_nuc('31P',diamond_tensor([7.585 2.942 2.328]*hz_per_mt,frame),[]);
    case 'np8'
        gmat=diamond_tensor([2.0016 2.0016 2.0048],diamond_frame_z([1 1 1]));
        nuclei{end+1}=diamond_nuc('31P',...
            diamond_tensor([3.2 3.2 5.6]*hz_per_mt,diamond_frame_z([1 1 1])),[]);
        nuclei{end+1}=diamond_nuc('31P',...
            diamond_tensor([8.8 8.8 13.6]*hz_per_mt,diamond_frame_z([1 1 1])),[]);
    case 'np9'
        gmat=diamond_tensor([2.0038 2.0038 2.0030],diamond_frame_z([1 1 1]));
        nuclei{end+1}=diamond_nuc('31P',...
            diamond_tensor([2.2 2.2 1.4]*hz_per_mt,diamond_frame_z([1 1 1])),[]);
    otherwise
        error('unknown phosphorus centre.');
end
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
if isfield(parameters,'include_13c')&&(~islogical(parameters.include_13c))
    error('parameters.include_13c must be logical.');
end
end

% The phosphorus table is compact, but each centre is separate data.

% Shared local helpers


% Field-unit conversion constants
function [hz_per_mt,hz_per_t]=diamond_hz_per_mt()
hz_per_mt=abs(spin('E'))/(2*pi)*1e-3;
hz_per_t=1e3*hz_per_mt;
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

