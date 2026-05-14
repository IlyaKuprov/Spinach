% Cobalt-related defect spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_co(parameters)
%
% Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),
% https://doi.org/10.3390/cryst7080237
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .centre       - 'o4' or 'nlo2'
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_co.m>

function [sys,inter]=diamond_co(parameters)

% Check input count
if nargin~=1
    error('exactly one input argument is required.');
end

% Check consistency
grumble(parameters);

% Set field-unit conversion constants
hz_per_mt=abs(spin('E'))/(2*pi)*1e-3;

% Select the cobalt centre
centre=lower(parameters.centre);
electron='E'; zfs=[];

% Get cobalt centre table data
switch centre
    case 'o4'
        gvals=[2.3463 1.8438 1.7045];
        Aco=[8.86 6.43 5.82];
        alpha=29;
    case 'nlo2'
        gvals=[2.3277 1.7982 1.7149];
        Aco=[8.24 6.57 5.76];
        alpha=28;
    otherwise
        error('unknown cobalt centre.');
end

% Assemble the cobalt principal-axis frame
yaxis=[0;1;1]/sqrt(2);
xbase=[1;0;0];
zbase=cross(xbase,yaxis);
xaxis=cosd(alpha)*xbase+sind(alpha)*zbase;
zaxis=cross(xaxis,yaxis);
frame=[xaxis(:) yaxis(:) zaxis(:)];

% Orthogonalise the frame
[frame,~]=qr(frame,0);

% Enforce right-handed axes
if det(frame)<0
    frame(:,3)=-frame(:,3);
end

% Build electron and cobalt tensors
gmat=((frame)*diag(gvals)*(frame)');
nuclei={struct('iso','59Co','A',((frame)*diag(Aco*hz_per_mt)*(frame)'),'Q',[])};

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
if ~isfield(parameters,'centre')
    error('parameters.centre field is required.');
end
if(~ischar(parameters.centre))
    error('parameters.centre must be a character string.');
end
if ~any(strcmpi(parameters.centre,{'o4','nlo2'}))
    error('parameters.centre must be ''o4'' or ''nlo2''.');
end
if ~isfield(parameters,'orientation')
    error('parameters.orientation field is required.');
end
if(~ischar(parameters.orientation))
    error('parameters.orientation must be a character string.');
end
if ~any(strcmp(parameters.orientation,{'111','110','100'}))
    error('parameters.orientation must be ''111'', ''110'', or ''100''.');
end
end

% No man is regular in his attendance at the House of
% Commons until he is married.
%
% Benjamin Disraeli

