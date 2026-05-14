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
%                      'r6', 'r10', or 'r11'
%
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_vacancy.m>

function [sys,inter]=diamond_vacancy(parameters)

% Check input count
if nargin~=1
    error('exactly one input argument is required.');
end

% Check consistency
grumble(parameters);

% Select the vacancy centre
centre=lower(parameters.centre);

% Set vacancy principal-axis frames
frame_111=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
            1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
            0          2/sqrt(6)  1/sqrt(3)];
frame_110=[-1/sqrt(2) 0 1/sqrt(2);...
            1/sqrt(2) 0 1/sqrt(2);...
            0         1 0];
switch centre
    case {'r4_w6','w6','r4'}
        electron='E3'; giso=2.0022;
        zfs=((frame_111)*diag([105 197 -303]*1e6)*(frame_111)');
    case 'w29'
        electron='E4'; giso=2.0019;
        zfs=((frame_111)*diag([297 156 -453]*1e6)*(frame_111)');
    case 'r5'
        electron='E3'; giso=2.0023;
        zfs=((frame_110)*diag([283 244 -524]*1e6)*(frame_110)');
    case 'o1'
        electron='E3'; giso=2.0023;
        zfs=((frame_110)*diag([109 95 -205]*1e6)*(frame_110)');
    case 'r6'
        electron='E3'; giso=2.0023;
        zfs=((frame_110)*diag([62 59 -120]*1e6)*(frame_110)');
    case 'r10'
        electron='E3'; giso=2.0023;
        zfs=((frame_110)*diag([36 36 -73]*1e6)*(frame_110)');
    case 'r11'
        electron='E3'; giso=2.0023;
        zfs=((frame_110)*diag([27 27 -53]*1e6)*(frame_110)');
    otherwise
        error('unknown vacancy centre.');
end

% Build the electron tensors
gmat=eye(3)*giso;
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
if ~isfield(parameters,'centre')
    error('parameters.centre must be specified.');
end
if ~isfield(parameters,'orientation')
    error('parameters.orientation must be specified.');
end
if ~ischar(parameters.centre)
    error('parameters.centre must be a character string.');
end
if ~ischar(parameters.orientation)
    error('parameters.orientation must be a character string.');
end
end


% A man is never so proud as when striking 
% an attitude of humility.
% 
% C.S. Lewis

