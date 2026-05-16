% NV0 excited-state spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_nv0_es(parameters)
%
% Magnetic parameters from Felton et al., Phys. Rev. B 77,
% 081201 (2008), https://doi.org/10.1103/PhysRevB.77.081201
%
% Parameters:
%
%    parameters is a structure with the following required fields:
%
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field
%      .nitrogen     - must be '15N'
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_nv0_es.m>

function [sys,inter]=diamond_nv0_es(parameters)

% Check input count
if nargin~=1
    error('exactly one input argument is required.');
end

% Check consistency
grumble(parameters);

% Build the electron tensors
frame=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
            1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
            0          2/sqrt(6)  1/sqrt(3)];
electron='E4';
gmat=((frame)*diag([2.0035 2.0035 2.0029])*(frame)');
zfs=frame*zfs2mat(1685e6,0,0,0,0)*frame';

% Add the nitrogen hyperfine tensor
nucleus=struct('iso','15N','A',((frame)*diag([-23.8e6 -23.8e6 -35.7e6])*(frame)'));

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
sys.isotopes={electron,nucleus.iso};
inter.zeeman.matrix=cell(1,numel(sys.isotopes));
inter.zeeman.matrix{1}=C*gmat*C';
inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes));
[~,~,zfs]=mat2ias(C*zfs*C');
inter.coupling.matrix{1,1}=zfs;
inter.coupling.matrix{1,2}=C*nucleus.A*C';

end

% Consistency enforcement
function grumble(parameters)
if(~isstruct(parameters))
    error('parameters must be a structure.');
end
if ~isfield(parameters,'orientation')
    error('parameters.orientation field is required.');
end
if(~ischar(parameters.orientation))
    error('parameters.orientation must be a character string.');
end
if ~ismember(parameters.orientation,{'111','110','100'})
    error('parameters.orientation must be ''111'', ''110'', or ''100''.');
end
if ~isfield(parameters,'nitrogen')
    error('parameters.nitrogen field is required.');
end
if(~ischar(parameters.nitrogen))
    error('parameters.nitrogen must be a character string.');
end
if ~strcmp(parameters.nitrogen,'15N')
    error('NV0 excited-state parameters are only available here for 15N.');
end
end

% In the fields of observation, chance 
% favors only the prepared mind.
%
% Louis Pasteur

