% SiV0 spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_siv0(parameters)
%
% Magnetic parameters from Edmonds et al., Phys. Rev. B 77,
% 245205 (2008), https://doi.org/10.1103/PhysRevB.77.245205
%
% Parameters:
%
%    parameters is a structure with the following required fields:
%
%      .silicon      - '29Si', 'none', or another silicon isotope
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field
%      .include_13c  - include reported 13C hyperfine couplings
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_siv0.m>

function [sys,inter]=diamond_siv0(parameters)

% Check input count
if nargin~=1
    error('exactly one input argument is required.');
end

% Check consistency
grumble(parameters);

% Build the electron tensors
electron='E3';
frame=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
            1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
            0          2/sqrt(6)  1/sqrt(3)];
gmat=((frame)*diag([2.0035 2.0035 2.0042])*(frame)');
zfs=frame*zfs2mat(1000e6,0,0,0,0)*frame';
nuclei={};

% Add the silicon isotope if requested
silicon=parameters.silicon;
if strcmp(silicon,'29Si')
    nuclei{end+1}=struct('iso','29Si','A',((frame)*diag([78.9e6 78.9e6 76.3e6])*(frame)'),'Q',[]);
elseif ~strcmp(silicon,'none')
    nuclei{end+1}=struct('iso',silicon,'A',zeros(3),'Q',[]);
end

% Add reported nearest-neighbour carbons
if parameters.include_13c
    Cmat=((frame)*diag([30.2e6 30.2e6 66.2e6])*(frame)');
    nuc_idx=numel(nuclei);
    nuclei(nuc_idx+1:nuc_idx+6)={[]};
    for n=1:6
        nuclei{nuc_idx+n}=struct('iso','13C','A',Cmat,'Q',[]);
    end
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
if ~isfield(parameters,'orientation')
    error('parameters.orientation field is missing.');
end
if ~ischar(parameters.orientation)
    error('parameters.orientation must be a character string.');
end
if ~ismember(parameters.orientation,{'111','110','100'})
    error('parameters.orientation must be ''111'', ''110'', or ''100''.');
end
if ~isfield(parameters,'silicon')
    error('parameters.silicon field is missing.');
end
if ~ischar(parameters.silicon)
    error('parameters.silicon must be a character string.');
end
if ~isfield(parameters,'include_13c')
    error('parameters.include_13c field is missing.');
end
if (~islogical(parameters.include_13c))||(~isscalar(parameters.include_13c))
    error('parameters.include_13c must be logical.');
end
end

% Among other evils which being unarmed brings you, 
% it causes you to be despised.
%
% Nicolo Machiavelli

