% Titanium-related defect spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_ti(parameters)
%
% Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),
% https://doi.org/10.3390/cryst7080237
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .centre       - 'n3' or 'ok1', default is 'n3'
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field, default is '111'
%      .titanium     - titanium isotope label, or 'none', default is '47Ti'
%      .include_13c - include reported 13C hyperfine couplings for OK1,
%                     false by default
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_ti.m>

function [sys,inter]=diamond_ti(parameters)

% Set default input
if nargin==0
    parameters=struct();
end

% Check consistency
grumble(parameters);

% Set default parameters
if ~isfield(parameters,'centre')
    parameters.centre='n3';
end
if ~isfield(parameters,'orientation')
    parameters.orientation='111';
end
if ~isfield(parameters,'include_13c')
    parameters.include_13c=false;
end

% Set field-unit conversion constants
hz_per_mt=abs(spin('E'))/(2*pi)*1e-3;

% Select the titanium centre
centre=lower(parameters.centre);
electron='E'; nuclei={}; zfs=[];
[gvals,An,Ati,g_alpha,A_alpha]=titanium_data(centre);
gframe=diamond_frame_alpha(g_alpha);
aframe=diamond_frame_alpha(A_alpha);
gmat=((gframe)*diag(gvals)*(gframe)');

% Add the nitrogen hyperfine tensor
nuclei{end+1}=struct('iso','14N','A',((aframe)*diag(An*hz_per_mt)*(aframe)'),'Q',[]);

% Add the titanium isotope if requested
if isfield(parameters,'titanium')
    titanium=parameters.titanium;
else
    titanium='47Ti';
end
if ~strcmp(titanium,'none')
    nuclei{end+1}=struct('iso',titanium,'A',((aframe)*diag(Ati*hz_per_mt)*(aframe)'),'Q',[]);
end

% Add reported OK1 nearest-neighbour carbons
if strcmp(centre,'ok1')&&parameters.include_13c
    Cmat=((diamond_frame_xz([1 1 0],[1 -1 -1]))*diag([2.62 2.62 4.38]*hz_per_mt)*(diamond_frame_xz([1 1 0],[1 -1 -1]))');
    nuclei{end+1}=struct('iso','13C','A',Cmat,'Q',[]);
    Cmat=((diamond_frame_xz([1 1 0],[-1 1 -1]))*diag([2.62 2.62 4.38]*hz_per_mt)*(diamond_frame_xz([1 1 0],[-1 1 -1]))');
    nuclei{end+1}=struct('iso','13C','A',Cmat,'Q',[]);
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

% Titanium centre table data
function [gvals,An,Ati,g_alpha,A_alpha]=titanium_data(centre)
switch centre
    case 'n3'
        gvals=[2.0022 2.0025 2.0020];
        An=[0.11 0.15 0.11];
        Ati=[0.28 0.40 0.28];
        g_alpha=32; A_alpha=26;
    case 'ok1'
        gvals=[2.0031 2.0019 2.0025];
        An=[0.55 0.77 0.54];
        Ati=[0.06 0.06 0.06];
        g_alpha=40; A_alpha=20;
    otherwise
        error('unknown titanium centre.');
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
if isfield(parameters,'titanium')&&(~ischar(parameters.titanium))
    error('parameters.titanium must be a character string.');
end
if isfield(parameters,'include_13c')&&(~islogical(parameters.include_13c))
    error('parameters.include_13c must be logical.');
end
end

% Make a principal-axis frame from three vectors
function frame=diamond_frame_xyz(xaxis,yaxis,zaxis)
frame=[xaxis(:) yaxis(:) zaxis(:)];
[frame,~]=qr(frame,0);
if det(frame)<0
    frame(:,3)=-frame(:,3);
end
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

% Principal-axis frame from a specified x axis and z axis
function frame=diamond_frame_xz(xaxis,zaxis)
xaxis=xaxis(:)/norm(xaxis,2);
zaxis=zaxis(:)-xaxis*dot(xaxis,zaxis(:));
zaxis=zaxis/norm(zaxis,2);
yaxis=cross(zaxis,xaxis);
frame=diamond_frame_xyz(xaxis,yaxis,zaxis);
end

