% Nickel-related defect spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_ni(parameters)
%
% W8 magnetic parameters from Ludwig and Woodbury, Phys. Rev. B 41,
% 3905 (1990), https://doi.org/10.1103/PhysRevB.41.3905
% Other nickel-centre table values from Nadolinny et al., Crystals
% 7, 237 (2017), https://doi.org/10.3390/cryst7080237
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .centre       - 'w8', 'ne1', 'ne2', 'ne3', 'ne4', 'ne5',
%                      'ne8', 'ab1', 'ab2', 'ab3', 'ab4', 'ab5',
%                      'nol1', or 'nirim5', default is 'w8'
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field, default is '111'
%      .nickel       - '61Ni', 'none', or another nickel isotope for W8,
%                      default is '61Ni'
%      .include_13c - include reported 13C hyperfine couplings,
%                     false by default
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_ni.m>

function [sys,inter]=diamond_ni(parameters)

% Set default input
if nargin==0
    parameters=struct();
end

% Check consistency
grumble(parameters);

% Set default parameters
if ~isfield(parameters,'centre')
    parameters.centre='w8';
end
if ~isfield(parameters,'orientation')
    parameters.orientation='111';
end
if ~isfield(parameters,'include_13c')
    parameters.include_13c=false;
end

% Set field-unit conversion constants
hz_per_mt=abs(spin('E'))/(2*pi)*1e-3;
hz_per_t=abs(spin('E'))/(2*pi);

% Select the nickel centre
centre=lower(parameters.centre);
nuclei={}; zfs=[];
switch centre
    case 'w8'
        electron='E4';
        gmat=eye(3)*2.032;
        if isfield(parameters,'nickel')
            nickel=parameters.nickel;
        else
            nickel='61Ni';
        end
        if strcmp(nickel,'61Ni')
            nuclei{end+1}=struct('iso','61Ni','A',eye(3)*0.65*hz_per_mt,'Q',[]);
        elseif ~strcmp(nickel,'none')
            nuclei{end+1}=struct('iso',nickel,'A',zeros(3),'Q',[]);
        end
        if parameters.include_13c
            Cmat=diamond_tensor([0.340 0.340 1.339]*hz_per_mt,diamond_frame_z([1 1 1]));
            nuc_idx=numel(nuclei);
            nuclei(nuc_idx+1:nuc_idx+4)={[]};
            for n=1:4
                nuclei{nuc_idx+n}=struct('iso','13C','A',Cmat,'Q',[]);
            end
        end
    case {'ne1','ne2','ne3','ne5','ne8'}
        electron='E';
        [gvals,avalues,alpha]=nickel_ne_data(centre);
        frame=diamond_frame_alpha(alpha);
        gmat=diamond_tensor(gvals,frame);
        nuclei=cell(1,size(avalues,1));
        for n=1:size(avalues,1)
            nuclei{n}=struct('iso','14N','A',diamond_tensor(avalues(n,:)*hz_per_mt,frame),'Q',[]);
        end
    case 'ne4'
        electron='E';
        frame=diamond_frame_z([1 1 1]);
        gmat=diamond_tensor([2.0988 2.0988 2.0227],frame);
    case {'ab1','ab2','ab3','ab4'}
        electron='E';
        [gvals,frame]=nickel_ab_data(centre);
        gmat=diamond_tensor(gvals,frame);
    case 'ab5'
        electron='E3';
        frame=diamond_frame_z([1 1 1]);
        gmat=diamond_tensor([2.022 2.022 2.037],frame);
        zfs=frame*zfs2mat(1.132*hz_per_t,0,0,0,0)*frame';
    case {'nol1','nirim5'}
        electron='E3';
        frame=diamond_frame_z([1 1 1]);
        gmat=diamond_tensor([2.002 2.002 2.0235],frame);
        zfs=frame*zfs2mat(-6.10*hz_per_t,0,0,0,0)*frame';
    otherwise
        error('unknown nickel centre.');
end

% Build the Spinach structures
[sys,inter]=diamond_system(electron,gmat,zfs,nuclei,parameters.orientation);

end

% Nickel-nitrogen centre table data
function [gvals,avalues,alpha]=nickel_ne_data(centre)
switch centre
    case 'ne1'
        gvals=[2.1282 2.0070 2.0908]; alpha=14;
        avalues=[2.09 1.43 1.45;2.09 1.43 1.45];
    case 'ne2'
        gvals=[2.1301 2.0100 2.0931]; alpha=14;
        avalues=[2.10 1.42 1.41;1.87 1.18 1.25;0.18 0.35 0.25];
    case 'ne3'
        gvals=[2.0729 2.0100 2.0476]; alpha=14;
        avalues=[1.60 1.24 1.15;0.66 0.50 0.50;0.66 0.50 0.50];
    case 'ne5'
        gvals=[2.0329 2.0898 2.0476]; alpha=27.5;
        avalues=[1.22 0.98 0.89;1.22 0.98 0.89];
    case 'ne8'
        gvals=[2.0439 2.1722 2.0476]; alpha=27.5;
        avalues=[1.14 0.78 0.75;1.14 0.78 0.75;1.14 0.78 0.75;1.14 0.78 0.75];
    otherwise
        error('unknown nickel NE centre.');
end
end

% Nitrogen-free nickel centre table data
function [gvals,frame]=nickel_ab_data(centre)
switch centre
    case 'ab1'
        frame=diamond_frame_z([1 1 1]); gvals=[2.0920 2.0920 2.0024];
    case 'ab2'
        frame=diamond_frame_z([1 1 1]); gvals=[2.0672 2.0672 2.0072];
    case 'ab3'
        frame=diamond_frame_xyz([1 0 0],[0 1 1],[0 -1 1]); gvals=[2.1105 2.0663 2.0181];
    case 'ab4'
        frame=diamond_frame_xyz([1 0 0],[0 1 1],[0 -1 1]); gvals=[2.0220 2.0094 2.0084];
    otherwise
        error('unknown nickel AB centre.');
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
if isfield(parameters,'nickel')&&(~ischar(parameters.nickel))
    error('parameters.nickel must be a character string.');
end
if isfield(parameters,'include_13c')&&(~islogical(parameters.include_13c))
    error('parameters.include_13c must be logical.');
end
end

% Nickel centres are a family, not one Hamiltonian.

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

% Make a principal-axis frame from three vectors
function frame=diamond_frame_xyz(xaxis,yaxis,zaxis)
frame=[xaxis(:) yaxis(:) zaxis(:)];
frame=diamond_frame_orth(frame);
end

% Frame used in the Nadolinny nickel and titanium tables
function frame=diamond_frame_alpha(alpha)
xaxis=[1;-1;0]/sqrt(2);
ybase=[1;1;0]/sqrt(2);
zbase=[0;0;1];
yaxis=cosd(alpha)*ybase+sind(alpha)*zbase;
zaxis=cross(xaxis,yaxis);
frame=diamond_frame_xyz(xaxis,yaxis,zaxis);
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

% Build the Spinach structures
function [sys,inter]=diamond_system(electron,gmat,zfs,nuclei,orientation)
C=diamond_orient(orientation);
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

