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
%                      'nol1', or 'nirim5'
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field
%      .nickel       - '61Ni', 'none', or another nickel isotope;
%                      required when .centre is 'w8'
%      .include_13c  - logical flag enabling reported 13C hyperfine couplings
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_ni.m>

function [sys,inter]=diamond_ni(parameters)

% Check input count
if nargin~=1
    error('exactly one input argument is required.');
end

% Check consistency
grumble(parameters);

% Set field-unit conversion constants
hz_per_mt=abs(spin('E'))/(2*pi)*1e-3;
hz_per_t=abs(spin('E'))/(2*pi);

% Select the nickel centre
centre=lower(parameters.centre);

% Set the trigonal principal-axis frame
frame_111=[-1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
            1/sqrt(2) -1/sqrt(6) 1/sqrt(3);...
            0          2/sqrt(6)  1/sqrt(3)];
nuclei={}; zfs=[];
switch centre
    case 'w8'
        electron='E4';
        gmat=eye(3)*2.032;
        nickel=parameters.nickel;
        if strcmp(nickel,'61Ni')
            nuclei{end+1}=struct('iso','61Ni','A',eye(3)*0.65*hz_per_mt);
        elseif ~strcmp(nickel,'none')
            nuclei{end+1}=struct('iso',nickel);
        end
        if parameters.include_13c
            Cmat=((frame_111)*diag([0.340 0.340 1.339]*hz_per_mt)*(frame_111)');
            nuc_idx=numel(nuclei);
            nuclei(nuc_idx+1:nuc_idx+4)={[]};
            for n=1:4
                nuclei{nuc_idx+n}=struct('iso','13C','A',Cmat);
            end
        end
    case {'ne1','ne2','ne3','ne5','ne8'}
        electron='E';

        % Get tabulated NE centre parameters
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
        end

        % Build the Nadolinny table frame
        xaxis=[1;-1;0]/sqrt(2);
        ybase=[1;1;0]/sqrt(2);
        zbase=[0;0;1];
        yaxis=cosd(alpha)*ybase+sind(alpha)*zbase;
        zaxis=cross(xaxis,yaxis);
        frame=diamond_frame_xyz(xaxis,yaxis,zaxis);
        gmat=((frame)*diag(gvals)*(frame)');
        nuclei=cell(1,size(avalues,1));
        for n=1:size(avalues,1)
            nuclei{n}=struct('iso','14N','A',((frame)*diag(avalues(n,:)*hz_per_mt)*(frame)'));
        end
    case 'ne4'
        electron='E';
        frame=frame_111;
        gmat=((frame)*diag([2.0988 2.0988 2.0227])*(frame)');
    case {'ab1','ab2','ab3','ab4'}
        electron='E';

        % Get tabulated AB centre parameters
        switch centre
            case 'ab1'
                frame=frame_111; gvals=[2.0920 2.0920 2.0024];
            case 'ab2'
                frame=frame_111; gvals=[2.0672 2.0672 2.0072];
            case 'ab3'
                frame=diamond_frame_xyz([1 0 0],[0 1 1],[0 -1 1]); gvals=[2.1105 2.0663 2.0181];
            case 'ab4'
                frame=diamond_frame_xyz([1 0 0],[0 1 1],[0 -1 1]); gvals=[2.0220 2.0094 2.0084];
        end
        gmat=((frame)*diag(gvals)*(frame)');
    case 'ab5'
        electron='E3';
        frame=frame_111;
        gmat=((frame)*diag([2.022 2.022 2.037])*(frame)');
        zfs=frame*zfs2mat(1.132*hz_per_t,0,0,0,0)*frame';
    case {'nol1','nirim5'}
        electron='E3';
        frame=frame_111;
        gmat=((frame)*diag([2.002 2.002 2.0235])*(frame)');
        zfs=frame*zfs2mat(-6.10*hz_per_t,0,0,0,0)*frame';
    otherwise
        error('unknown nickel centre.');
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
if ~isempty(zfs)
    [~,~,zfs]=mat2ias(C*zfs*C');
    inter.coupling.matrix{1,1}=zfs;
end
for n=1:numel(nuclei)
    if isfield(nuclei{n},'A')&&norm(nuclei{n}.A,2)>0
        inter.coupling.matrix{1,n+1}=C*nuclei{n}.A*C';
    end
end

end

% Consistency enforcement
function grumble(parameters)
if ~isstruct(parameters)
    error('parameters must be a structure.');
end
required={'centre','orientation','include_13c'};
for n=1:numel(required)
    if ~isfield(parameters,required{n})
        error(['parameters.' required{n} ' field is required.']);
    end
end
if ~ischar(parameters.centre)
    error('parameters.centre must be a character string.');
end
if ~ischar(parameters.orientation)
    error('parameters.orientation must be a character string.');
end
if ~islogical(parameters.include_13c)||~isscalar(parameters.include_13c)
    error('parameters.include_13c must be a scalar logical.');
end
if strcmpi(parameters.centre,'w8')&&(~isfield(parameters,'nickel'))
    error('parameters.nickel field is required for W8 centre.');
end
if isfield(parameters,'nickel')&&(~ischar(parameters.nickel))
    error('parameters.nickel must be a character string.');
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

% I have brought myself, by long meditation, to the conviction 
% that a human being with a settled purpose must accomplish it, 
% and that nothing can resist a will which would stake even 
% existence upon its fulfillment. 
%
% Benjamin Disraeli

