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
%                      'np6', 'np8', or 'np9'
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field
%      .include_13c - include reported 13C hyperfine couplings for MA1
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_p.m>

function [sys,inter]=diamond_p(parameters)

% Check input count
if nargin~=1
    error('exactly one input argument is required.');
end

% Check consistency
grumble(parameters);

% Set field-unit conversion constants
hz_per_mt=abs(spin('E'))/(2*pi)*1e-3;

% Select the phosphorus centre
centre=lower(parameters.centre);
electron='E';
nuclei={}; frame=eye(3);

% Set the trigonal principal-axis frame
frame_111=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
            1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
            0          2/sqrt(6)  1/sqrt(3)];

% Get phosphorus centre table data
switch centre
    case 'ma1'
        gmat=eye(3)*2.0025;
        nuclei{end+1}=struct('iso','31P','A',...
            ((frame_111)*diag([1.96 1.96 2.32]*hz_per_mt)*(frame_111)'));
        if parameters.include_13c
            nuclei{end+1}=struct('iso','13C','A',...
                ((frame_111)*diag([13.92 13.92 18.13]*hz_per_mt)*(frame_111)'));
        end
    case 'np1'
        gmat=((frame)*diag([2.00243 2.0028 2.0026])*(frame)');
        nuclei{end+1}=struct('iso','31P','A',((frame)*diag([2.08 2.02 2.18]*hz_per_mt)*(frame)'));
        nuclei{end+1}=struct('iso','14N','A',((frame)*diag([4.08 3.10 3.00]*hz_per_mt)*(frame)'));
    case 'np2'
        gmat=eye(3)*2.0025;
        nuclei{end+1}=struct('iso','31P','A',...
            ((frame_111)*diag([2.09 2.09 2.34]*hz_per_mt)*(frame_111)'));
        nuclei{end+1}=struct('iso','14N','A',...
            ((frame_111)*diag([3.09 3.09 6.42]*hz_per_mt)*(frame_111)'));
    case 'np3'
        gmat=eye(3)*2.0025;
        nuclei{end+1}=struct('iso','31P','A',...
            ((frame_111)*diag([18.23 18.23 17.48]*hz_per_mt)*(frame_111)'));
        nuclei{end+1}=struct('iso','14N','A',...
            ((frame_111)*diag([0.33 0.33 0.10]*hz_per_mt)*(frame_111)'));
    case 'np4'
        gmat=((frame)*diag([2.0009 2.0012 2.00047])*(frame)');
        nuclei{end+1}=struct('iso','31P','A',((frame)*diag([5.456 3.838 3.80]*hz_per_mt)*(frame)'));
    case 'np5'
        gmat=((frame_111)*diag([2.0009 2.0009 2.00087])*(frame_111)');
        nuclei{end+1}=struct('iso','31P','A',...
            ((frame_111)*diag([1.024 1.024 6.522]*hz_per_mt)*(frame_111)'));
    case 'np6'
        gmat=((frame_111)*diag([2.00083 2.00083 2.00085])*(frame_111)');
        nuclei{end+1}=struct('iso','31P','A',((frame)*diag([7.585 2.942 2.328]*hz_per_mt)*(frame)'));
    case 'np8'
        gmat=((frame_111)*diag([2.0016 2.0016 2.0048])*(frame_111)');
        nuclei{end+1}=struct('iso','31P','A',...
            ((frame_111)*diag([3.2 3.2 5.6]*hz_per_mt)*(frame_111)'));
        nuclei{end+1}=struct('iso','31P','A',...
            ((frame_111)*diag([8.8 8.8 13.6]*hz_per_mt)*(frame_111)'));
    case 'np9'
        gmat=((frame_111)*diag([2.0038 2.0038 2.0030])*(frame_111)');
        nuclei{end+1}=struct('iso','31P','A',...
            ((frame_111)*diag([2.2 2.2 1.4]*hz_per_mt)*(frame_111)'));
    otherwise
        error('unknown phosphorus centre.');
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
if ~isfield(parameters,'centre')
    error('parameters.centre field is required.');
end
if(~ischar(parameters.centre))
    error('parameters.centre must be a character string.');
end
if ~ismember(lower(parameters.centre),{'ma1','np1','np2','np3','np4','np5','np6','np8','np9'})
    error('parameters.centre has an unsupported value.');
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
if ~isfield(parameters,'include_13c')
    error('parameters.include_13c field is required.');
end
if(~islogical(parameters.include_13c))
    error('parameters.include_13c must be logical.');
end
end

% The fact that we live at the bottom of a deep gravity well, on the
% surface of a gas covered planet going around a nuclear fireball 90
% million miles away and think this to be normal is obviously some
% indication of how skewed our perspective tends to be.
%
% Douglas Adams

