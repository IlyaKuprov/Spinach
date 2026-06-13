% Nitrogen interstitial spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_n_inter(parameters)
%
% Magnetic parameters from Felton et al., J. Phys. Condens. Matter
% 21, 364212 (2009), https://doi.org/10.1088/0953-8984/21/36/364212
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .centre       - 'war9' or 'war10'
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field
%      .nitrogen     - '14N' or '15N'
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% alexey.bogdanov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=diamond_n_inter.m>

function [sys,inter]=diamond_n_inter(parameters)

% Check input count
if nargin~=1
    error('exactly one input argument is required.');
end

% Check consistency
grumble(parameters);

% Build the common frame
electron='E';
gframe=diamond_frame_xyz([sind(90)*cosd(45);sind(90)*sind(45);cosd(90)],...
                         [sind(180)*cosd(45);sind(180)*sind(45);cosd(180)],...
                         [sind(90)*cosd(315);sind(90)*sind(315);cosd(90)]);

% Select the interstitial centre
centre=lower(parameters.centre);
switch centre
    case 'war9'
        gmat=((gframe)*diag([2.00343 2.00272 2.00268])*(gframe)');
        Amat=((gframe)*diag([8.30e6 7.85e6 8.17e6])*(gframe)');
    case 'war10'
        aframe=diamond_frame_xyz([sind(44.8)*cosd(45.0);sind(44.8)*sind(45.0);cosd(44.8)],...
                                 [sind(134.8)*cosd(45.0);sind(134.8)*sind(45.0);cosd(134.8)],...
                                 [sind(90)*cosd(315);sind(90)*sind(315);cosd(90)]);
        gmat=((gframe)*diag([2.00344 2.00272 2.00269])*(gframe)');
        Amat=((aframe)*diag([1.00e6 -1.01e6 0.00e6])*(aframe)');
    otherwise
        error('unknown nitrogen interstitial centre.');
end

% Add the nitrogen hyperfine tensor
switch parameters.nitrogen
    case '15N'
        nucleus=struct('iso','15N','A',Amat);
    case '14N'
        scale=spin('14N')/spin('15N');
        nucleus=struct('iso','14N','A',scale*Amat);
    otherwise
        error('parameters.nitrogen must be ''14N'' or ''15N''.');
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
sys.isotopes={electron,nucleus.iso};
inter.zeeman.matrix=cell(1,numel(sys.isotopes));
inter.zeeman.matrix{1}=C*gmat*C';
inter.coupling.matrix=cell(numel(sys.isotopes),numel(sys.isotopes));
inter.coupling.matrix{1,2}=C*nucleus.A*C';

end

% Make a principal-axis frame from three vectors
function frame=diamond_frame_xyz(xaxis,yaxis,zaxis)
frame=[xaxis(:) yaxis(:) zaxis(:)];
[frame,~]=qr(frame,0);
if det(frame)<0
    frame(:,3)=-frame(:,3);
end
end

% Consistency enforcement
function grumble(parameters)
if(~isstruct(parameters))
    error('parameters must be a structure.');
end
if(~isfield(parameters,'centre'))
    error('parameters.centre field is required.');
end
if(~ischar(parameters.centre))
    error('parameters.centre must be a character string.');
end
if(~isfield(parameters,'orientation'))
    error('parameters.orientation field is required.');
end
if(~ischar(parameters.orientation))
    error('parameters.orientation must be a character string.');
end
if(~isfield(parameters,'nitrogen'))
    error('parameters.nitrogen field is required.');
end
if(~ischar(parameters.nitrogen))
    error('parameters.nitrogen must be a character string.');
end
if(~any(strcmp(parameters.nitrogen,{'14N','15N'})))
    error('parameters.nitrogen must be ''14N'' or ''15N''.');
end
end

% It is a general rule of human nature that people despise
% those who treat them well, and look up to those who make
% no concessions. 
%
% Thucydides

