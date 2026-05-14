% N2V- spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_n2vm(parameters)
%
% Magnetic parameters from Green et al., Phys. Rev. B 92,
% 165204 (2015), https://doi.org/10.1103/PhysRevB.92.165204
%
% Parameters:
%
%    parameters is a structure with the following required fields:
%
%      .nitrogen     - '14N' or '15N'
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field
%      .include_13c  - include reported 13C hyperfine couplings,
%                      true or false
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_n2vm.m>

function [sys,inter]=diamond_n2vm(parameters)

% Check input count
if nargin~=1
    error('exactly one input argument is required.');
end

% Check consistency
grumble(parameters);

% Build the reported frames
c2rot=anax2dcm([0 0 1],-pi);
gframe=diamond_frame_xyz([1 1 0],[0 0 1],[1 -1 0]);
nframe=anax2dcm([1 -1 0],3.5*pi/180)*...
       diamond_frame_xyz([1 1 -2],[1 1 1],[1 -1 0]);

% Orthogonalise the carbon hyperfine frame axes
xaxis=[-1 -1 0]';
xaxis=xaxis/norm(xaxis,2);
zaxis=[-1 1 1]';
zaxis=zaxis-xaxis*dot(xaxis,zaxis);
zaxis=zaxis/norm(zaxis,2);
yaxis=cross(zaxis,xaxis);
cframe=anax2dcm([-1 -1 0],-2.0*pi/180)*...
       diamond_frame_xyz(xaxis,yaxis,zaxis);

% Build the electron tensors
electron='E';
gmat=((gframe)*diag([2.00345 2.00274 2.00271])*(gframe)');
zfs=[]; nuclei={};

% Add nitrogen tensors
if strcmp(parameters.nitrogen,'15N')
    nuclei{end+1}=struct('iso','15N','A',((nframe)*diag([3.47e6 4.51e6 4.09e6])*(nframe)'),'Q',[]);
    nuclei{end+1}=struct('iso','15N','A',((c2rot*nframe)*diag([3.47e6 4.51e6 4.09e6])*(c2rot*nframe)'),'Q',[]);
elseif strcmp(parameters.nitrogen,'14N')
    scale=spin('14N')/spin('15N');
    Qframe=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
            1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
            0          2/sqrt(6)  1/sqrt(3)];
    nuclei{end+1}=struct('iso','14N','A',scale*((nframe)*diag([3.47e6 4.51e6 4.09e6])*(nframe)'),'Q',...
                               Qframe*zfs2mat(-5.0e6,0,0,0,0)*Qframe');
    nuclei{end+1}=struct('iso','14N','A',scale*((c2rot*nframe)*diag([3.47e6 4.51e6 4.09e6])*(c2rot*nframe)'),'Q',...
                               ((c2rot*Qframe)*zfs2mat(-5.0e6,0,0,0,0)*(c2rot*Qframe)'));
else
    error('parameters.nitrogen must be ''14N'' or ''15N''.');
end

% Add reported nearest-neighbour carbons
if parameters.include_13c
    nuclei{end+1}=struct('iso','13C','A',((cframe)*diag([202.3e6 202.3e6 317.5e6])*(cframe)'),'Q',[]);
    nuclei{end+1}=struct('iso','13C','A',((c2rot*cframe)*diag([202.3e6 202.3e6 317.5e6])*(c2rot*cframe)'),'Q',[]);
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

% Consistency enforcement
function grumble(parameters)
if(~isstruct(parameters))
    error('parameters must be a structure.');
end
if(~isfield(parameters,'orientation'))
    error('parameters.orientation field is required.');
end
if(~isfield(parameters,'nitrogen'))
    error('parameters.nitrogen field is required.');
end
if(~isfield(parameters,'include_13c'))
    error('parameters.include_13c field is required.');
end
if(~ischar(parameters.orientation))
    error('parameters.orientation must be a character string.');
end
if(~any(strcmp(parameters.orientation,{'111','110','100'})))
    error('parameters.orientation must be ''111'', ''110'', or ''100''.');
end
if(~ischar(parameters.nitrogen))
    error('parameters.nitrogen must be a character string.');
end
if(~any(strcmp(parameters.nitrogen,{'14N','15N'})))
    error('parameters.nitrogen must be ''14N'' or ''15N''.');
end
if(~islogical(parameters.include_13c)||~isscalar(parameters.include_13c))
    error('parameters.include_13c must be a logical scalar.');
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

% If nothing is self-evident, nothing can be proved.
% 
% C.S. Lewis

