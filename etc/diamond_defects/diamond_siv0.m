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
%      .n_13c       - number of reported nearest-neighbour 13C
%                     hyperfine couplings, between 0 and 6
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
    nuclei{end+1}=struct('iso','29Si','A',((frame)*diag([78.9e6 78.9e6 76.3e6])*(frame)'));
elseif ~strcmp(silicon,'none')
    nuclei{end+1}=struct('iso',silicon);
end

% Add reported nearest-neighbour carbons
if parameters.n_13c>0
    Cmat=((frame)*diag([30.2e6 30.2e6 66.2e6])*(frame)');
    nuc_idx=numel(nuclei);
    nuclei(nuc_idx+1:nuc_idx+parameters.n_13c)={[]};
    for n=1:parameters.n_13c
        nuclei{nuc_idx+n}=struct('iso','13C','A',Cmat);
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
for n=1:numel(nuclei)
    sys.isotopes{n+1}=nuclei{n}.iso;
end
inter.zeeman.matrix=cell(1,numel(sys.isotopes));
inter.zeeman.matrix{1}=C*gmat*C';
inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes));
[~,~,zfs]=mat2ias(C*zfs*C');
inter.coupling.matrix{1,1}=zfs;
for n=1:numel(nuclei)
    if isfield(nuclei{n},'A')&&norm(nuclei{n}.A,2)>0
        inter.coupling.matrix{1,n+1}=C*nuclei{n}.A*C';
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
if ~isfield(parameters,'n_13c')
    error('parameters.n_13c field is missing.');
end
if (~isnumeric(parameters.n_13c))||(~isscalar(parameters.n_13c))||...
   (parameters.n_13c<0)||(parameters.n_13c>6)||(mod(parameters.n_13c,1)~=0)
    error('parameters.n_13c must be an integer between 0 and 6.');
end
end

% Among other evils which being unarmed brings you, 
% it causes you to be despised.
%
% Nicolo Machiavelli

