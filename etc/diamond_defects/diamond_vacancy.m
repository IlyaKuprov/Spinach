% Vacancy-family defect spin systems for diamond. Syntax:
%
%          [sys,inter]=diamond_vacancy(parameters)
%
% R4/W6 parameters from Twitchen et al., Phys. Rev. B 59, 12900
% (1999), https://doi.org/10.1103/PhysRevB.59.12900; W29
% parameters from Kirui et al., Diam. Relat. Mater. 8, 1569
% (1999), https://doi.org/10.1016/S0925-9635(99)00037-0;
% R5/O1/R6/R10/R11 parameters from Iakoubovskii and Stesmans,
% Phys. Rev. B 66, 045406 (2002),
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

% Set R4/W6 principal-axis frame
r4_x=[1;-1;0]/sqrt(2);
r4_z=[sind(54.2)*cosd(45);sind(54.2)*sind(45);cosd(54.2)];
r4_z=r4_z/norm(r4_z);
r4_y=cross(r4_z,r4_x);
r4_y=r4_y/norm(r4_y);
frame_r4=[r4_x r4_y r4_z];

% Set W29 principal-axis frame
w29_x=[0;-1;1]/sqrt(2);
w29_z=[0.619;-0.556;-0.556];
w29_z=w29_z/norm(w29_z);
w29_y=cross(w29_z,w29_x);
w29_y=w29_y/norm(w29_y);
frame_w29=[w29_x w29_y w29_z];

% Set vacancy-chain principal-axis frame
frame_110=[-1/sqrt(2) 0 1/sqrt(2);...
            1/sqrt(2) 0 1/sqrt(2);...
            0         1 0];
switch centre
    case {'r4_w6','w6','r4'}
        electron='E3';
        gmat=((frame_r4)*diag([2.0022 2.0026 2.0013])*(frame_r4)');
        zfs=((frame_r4)*diag([105 197 -303]*1e6)*(frame_r4)');
    case 'w29'
        electron='E4';
        gmat=((frame_w29)*diag([2.002 1.997 2.005])*(frame_w29)');
        zfs=((frame_w29)*diag([297 156 -453]*1e6)*(frame_w29)');
    case 'r5'
        electron='E3';
        gmat=((frame_110)*diag([2.00275 2.00265 2.00205])*(frame_110)');
        zfs=((frame_110)*diag([283 244 -524]*1e6)*(frame_110)');
    case 'o1'
        electron='E3';
        gmat=((frame_110)*diag([2.00299 2.00273 2.00212])*(frame_110)');
        zfs=((frame_110)*diag([109 95 -205]*1e6)*(frame_110)');
    case 'r6'
        electron='E3';
        gmat=((frame_110)*diag([2.00299 2.00273 2.00212])*(frame_110)');
        zfs=((frame_110)*diag([62 59 -120]*1e6)*(frame_110)');
    case 'r10'
        electron='E3';
        gmat=((frame_110)*diag([2.00295 2.00269 2.00212])*(frame_110)');
        zfs=((frame_110)*diag([36 36 -73]*1e6)*(frame_110)');
    case 'r11'
        electron='E3';
        gmat=((frame_110)*diag([2.00301 2.00278 2.00])*(frame_110)');
        zfs=((frame_110)*diag([27 27 -53]*1e6)*(frame_110)');
    otherwise
        error('unknown vacancy centre.');
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
inter.zeeman.matrix=cell(1,numel(sys.isotopes));
inter.zeeman.matrix{1}=C*gmat*C';
inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes));
[~,~,zfs]=mat2ias(C*zfs*C');
inter.coupling.matrix{1,1}=zfs;

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
