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
%      .centre       - 'war9' or 'war10', default is 'war9'
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field, default is '111'
%      .nitrogen     - must be '15N', default is '15N'
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_n_inter.m>

function [sys,inter]=diamond_n_inter(parameters)

% Set default input
if nargin==0
    parameters=struct();
end

% Check consistency
grumble(parameters);

% Set default parameters
if ~isfield(parameters,'centre')
    parameters.centre='war9';
end
if ~isfield(parameters,'orientation')
    parameters.orientation='111';
end
if ~isfield(parameters,'nitrogen')
    parameters.nitrogen='15N';
end

% Restrict to the verified isotope
if ~strcmp(parameters.nitrogen,'15N')
    error('WAR9/WAR10 parameters are only available here for 15N.');
end

% Build the common frame
electron='E'; zfs=[]; nuclei={};
gframe=diamond_frame_xyz(diamond_sph_vec(90,45),diamond_sph_vec(180,45),diamond_sph_vec(90,315));

% Select the interstitial centre
centre=lower(parameters.centre);
switch centre
    case 'war9'
        gmat=diamond_tensor([2.00343 2.00272 2.00268],gframe);
        Amat=diamond_tensor([8.30e6 7.85e6 8.17e6],gframe);
    case 'war10'
        aframe=diamond_frame_xyz(diamond_sph_vec(44.8,45.0),...
                                 diamond_sph_vec(134.8,45.0),...
                                 diamond_sph_vec(90,315));
        gmat=diamond_tensor([2.00344 2.00272 2.00269],gframe);
        Amat=diamond_tensor([1.00e6 -1.01e6 0.00e6],aframe);
    otherwise
        error('unknown nitrogen interstitial centre.');
end

% Add the nitrogen hyperfine tensor
nuclei{end+1}=diamond_nuc('15N',Amat,[]);

% Build the Spinach structures
[sys,inter]=diamond_system(electron,gmat,zfs,nuclei,parameters.orientation);

end

% Consistency enforcement
function grumble(parameters)
if(~isstruct(parameters))
    error('parameters must be a structure.');
end
if isfield(parameters,'centre')&&(~ischar(parameters.centre))
    error('parameters.centre must be a character string.');
end
if isfield(parameters,'orientation')&&(~ischar(parameters.orientation))
    error('parameters.orientation must be a character string.');
end
if isfield(parameters,'nitrogen')&&(~ischar(parameters.nitrogen))
    error('parameters.nitrogen must be a character string.');
end
end

% Interstitial names should not hide their isotope assumptions.

