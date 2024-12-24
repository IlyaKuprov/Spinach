% Crops 2D spectra to user-specified ranges (in ppm), respecting the
% digital resolution. Syntax:
%
% [spec,parameters]=crop_2d(spin_system,spec,parameters,crop_ranges)
%
% Inputs:
%
%          spec - 2D matrix containing the spectrum
%
%   crop_ranges - cropping bounds, supplied in the following
%                 format: {[f1_min f1_max],[f2_min f2_max]}
%
% The following subfields are required in the parameters structure:
%
%    parameters.sweep   -    one or two sweep widths, Hz
%
%    parameters.spins   -    cell array with one ot two character
%                            strings specifying the working spins.
%
%    parameters.offset  -    one or two transmitter offsets, Hz
%
% Outputs:
%
%          spec - 2D matrix containing the cropped spectrum
%
%    parameters - the updated parameters structure
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=crop_2d.m>

function [spec,parameters]=crop_2d(spin_system,spec,parameters,crop_ranges)

% Check consistency
grumble(spin_system,spec,parameters,crop_ranges);

% Accommodate homonuclear 2D sequences
if isscalar(parameters.spins)
    parameters.spins= [parameters.spins  parameters.spins];
end
if isscalar(parameters.offset)
    parameters.offset=[parameters.offset parameters.offset];
end
if isscalar(parameters.sweep)
    parameters.sweep= [parameters.sweep  parameters.sweep];
end

% Build axes and apply offsets
axis_f1_hz=linspace(-parameters.sweep(1)/2,parameters.sweep(1)/2,size(spec,1))+parameters.offset(1);
axis_f2_hz=linspace(-parameters.sweep(2)/2,parameters.sweep(2)/2,size(spec,2))+parameters.offset(2);

% Convert the units 
axis_f1_ppm=1e6*(2*pi)*axis_f1_hz/(spin(parameters.spins{1})*spin_system.inter.magnet);
axis_f2_ppm=1e6*(2*pi)*axis_f2_hz/(spin(parameters.spins{2})*spin_system.inter.magnet);

% Find array bounds
l_bound_f1=find(axis_f1_ppm>crop_ranges{1}(1),1);
r_bound_f1=find(axis_f1_ppm>crop_ranges{1}(2),1);
l_bound_f2=find(axis_f2_ppm>crop_ranges{2}(1),1);
r_bound_f2=find(axis_f2_ppm>crop_ranges{2}(2),1);

% Find the new offsets
parameters.offset=[(axis_f1_hz(l_bound_f1)+axis_f1_hz(r_bound_f1))/2 ...
                   (axis_f2_hz(l_bound_f2)+axis_f2_hz(r_bound_f2))/2];

% Find the new sweeps
parameters.sweep= [(axis_f1_hz(r_bound_f1)-axis_f1_hz(l_bound_f1)) ...
                   (axis_f2_hz(r_bound_f2)-axis_f2_hz(l_bound_f2))];
               
% Update the point counts
parameters.zerofill=[(r_bound_f1-l_bound_f1) (r_bound_f2-l_bound_f2)];

% Cut the spectrum
spec=spec(l_bound_f1:r_bound_f1,l_bound_f2:r_bound_f2);

end

% Consistency enforcement
function grumble(spin_system,spec,parameters,crop_ranges)
if (~isnumeric(spec))||(~ismatrix(spec))
    error('spec must be a matrix.');
end
if (~isfield(parameters,'offset'))
    error('offsets should be specified in parameters.offset variable.');
end
if (numel(parameters.offset)~=1)&&(numel(parameters.offset)~=2)
    error('parameters.offset array should have one or two elements.');
end
if ~isfield(parameters,'sweep')
    error('sweep widths should be specified in parameters.sweep variable.');
end
if (numel(parameters.sweep)~=1)&&(numel(parameters.sweep)~=2)
    error('parameters.sweep array should have one or two elements.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
end
if ~iscell(parameters.spins)
    error('parameters.spins should be a cell array of character strings.');
end
if (numel(parameters.spins)~=1)&&(numel(parameters.spins)~=2)
    error('parameters.spins cell array should have one or two elements.');
end
for n=1:numel(parameters.spins)
    if ~ismember(parameters.spins{n},spin_system.comp.isotopes)
        error('one of the isotopes in parameters.spins is not present in the system.');
    end
end
if (~iscell(crop_ranges))||(numel(crop_ranges)~=2)
    error('crop_ranges must be a cell array with two vectors.');
end
if (crop_ranges{1}(1)>=crop_ranges{1}(2))||(crop_ranges{2}(1)>=crop_ranges{2}(2))
    error('the vectors in the crop_ranges array must have their elements in ascending order.');
end
end

% There is long term price instability for lanthanides, as most come 
% from China, but are produced with little regard for the environment
% or workers' health. Attempts to improve conditions led to a 3000%
% price increase for dysprosium in 2011.
%
% EPSRC reviewer, on one of IK's grant applications

