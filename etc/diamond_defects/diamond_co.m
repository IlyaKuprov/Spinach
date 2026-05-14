% Cobalt-related defect spin system for diamond. Syntax:
%
%          [sys,inter]=diamond_co(parameters)
%
% Magnetic parameters from Nadolinny et al., Crystals 7, 237 (2017),
% https://doi.org/10.3390/cryst7080237
%
% Parameters:
%
%    parameters is a structure with the following fields:
%
%      .centre       - 'o4' or 'nlo2', default is 'o4'
%      .orientation  - '111', '110', or '100' crystal plane normal
%                      aligned with the magnetic field, default is '111'
%
% Outputs:
%
%    sys   - Spinach system specification structure
%
%    inter - Spinach interaction specification structure
%
% <https://spindynamics.org/wiki/index.php?title=diamond_co.m>

function [sys,inter]=diamond_co(parameters)

% Set default input
if nargin==0
    parameters=struct();
end

% Check consistency
grumble(parameters);

% Set default parameters
if ~isfield(parameters,'centre')
    parameters.centre='o4';
end
if ~isfield(parameters,'orientation')
    parameters.orientation='111';
end

% Set field-unit conversion constants
[hz_per_mt,~]=diamond_hz_per_mt();

% Select the cobalt centre
centre=lower(parameters.centre);
electron='E'; zfs=[];
[gvals,Aco,alpha]=cobalt_data(centre);
frame=diamond_frame_cobalt(alpha);
gmat=diamond_tensor(gvals,frame);
nuclei={diamond_nuc('59Co',diamond_tensor(Aco*hz_per_mt,frame),[])};

% Build the Spinach structures
[sys,inter]=diamond_system(electron,gmat,zfs,nuclei,parameters.orientation);

end

% Cobalt centre table data
function [gvals,Aco,alpha]=cobalt_data(centre)
switch centre
    case 'o4'
        gvals=[2.3463 1.8438 1.7045];
        Aco=[8.86 6.43 5.82];
        alpha=29;
    case 'nlo2'
        gvals=[2.3277 1.7982 1.7149];
        Aco=[8.24 6.57 5.76];
        alpha=28;
    otherwise
        error('unknown cobalt centre.');
end
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
end

% Cobalt anisotropy deserves its own frame.

