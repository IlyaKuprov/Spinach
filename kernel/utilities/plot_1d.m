% 1D plotting utility. Syntax:
%
%           plot_1d(spin_system,spectrum,parameters,varargin)
%
% Parameters:
%
%    spectrum                      a column vector containing the
%                                  spectrum
%
%    parameters.sweep              sweep width, Hz
%
%    parameters.spins              spin species, e.g. {'1H'}
%
%    parameters.offset             transmitter offset, Hz
%
%    parameters.axis_units         axis units ('ppm','Gauss',
%                                  'mT','T','Hz','kHz','MHz')
%
%    parameters.derivative         if set to 1, the spectrum is
%                                  differentiated before plotting
%
%    parameters.invert_axis        if set to 1, the frequency axis 
%                                  is inverted before plotting
%
%    varargin                      any number of any other para-
%                                  meters; these will be passed
%                                  to Matlab's plot() function
%
% Outputs:
%
%    this function produces a figure
%
% i.kuprov@soton.ac.uk
% ledwards@cbs.mpg.de
%
% <https://spindynamics.org/wiki/index.php?title=plot_1d.m>

function plot_1d(spin_system,spectrum,parameters,varargin)

% Set common defaults
parameters=defaults(spin_system,parameters);

% Check consistency
grumble(spin_system,spectrum,parameters);

% If a complex spectrum is received, plot both components
if ~isreal(spectrum)
    
    % Recursively plot the real and imaginary components
    plot_1d(spin_system,real(spectrum),parameters,varargin{:}); hold on;
    plot_1d(spin_system,imag(spectrum),parameters,varargin{:});
    legend({'real','imag'},'Location','NorthEast'); return;
    
end

% Get the axis
[ax,ax_label]=axis_1d(spin_system,parameters);

% Compute the derivative if necessary
if isfield(parameters,'derivative')&&parameters.derivative
    spectrum=fdvec(spectrum,5,1);
end

% Plot the spectrum
plot(ax,spectrum,varargin{:}); 
xlim tight; ylim padded;
box on; kgrid;

% Label the axis
kxlabel(ax_label);

% Invert the axis if necessary
if isfield(parameters,'invert_axis')&&parameters.invert_axis
    set(gca,'XDir','reverse');
end

end

% Default parameters
function parameters=defaults(spin_system,parameters)
if (~isfield(parameters,'offset'))&&isscalar(parameters.sweep)
    report(spin_system,'parameters.offset field not set, assuming zero offsets.');
    parameters.offset=zeros(size(parameters.spins));
end
if ~isfield(parameters,'axis_units')
    report(spin_system,'parameters.axis_units field not set, assuming ppm.');
    parameters.axis_units='ppm';
end
if ~isfield(parameters,'invert_axis')
    report(spin_system,'parameters.invert_axis field not set, assuming NMR tradition.');
    parameters.invert_axis=1;
end
if (~isfield(parameters,'zerofill'))&&isfield(parameters,'npoints')
    parameters.zerofill=parameters.npoints;
end
end

% Consistency enforcement
function grumble(spin_system,spectrum,parameters)
if (~isnumeric(spectrum))
    error('spectrum must be a numeric array.');
end
if ~isfield(parameters,'sweep')
    error('sweep width should be specified in parameters.sweep variable.');
end
if (~isnumeric(parameters.sweep))||(~isreal(parameters.sweep))
    error('elements of parameters.sweep must be real numbers.');
end
if numel(parameters.sweep)>2
    error('parameters.sweep should have one or two elements.');
end
if numel(parameters.sweep)==2
    if (parameters.sweep(1)>=parameters.sweep(2))
        error('the first element of parameters.sweep must be smaller than the second.'); 
    end
    if isfield(parameters,'offset')
        error('offset should not be specified when parameters.sweep has two elements.');
    end
else
    if ~isfield(parameters,'offset')
        error('offset should be specified in parameters.offset variable.');
    end
    if numel(parameters.offset)~=1
        error('parameters.offset vector should have exactly one element.');
    end
end
if ~isfield(parameters,'axis_units')
    error('axis units must be specified in parameters.axis_units variable.');
end
if ~ischar(parameters.axis_units)
    error('parameters.axis_units must be a character string.');
end
if ~isfield(parameters,'spins')
    error('working spins should be specified in parameters.spins variable.');
end
if (~iscell(parameters.spins))||(numel(parameters.spins)~=1)
    error('parameters.spins must be a cell array with exactly one element.');
end
if ~ischar(parameters.spins{1})
    error('elements of parameters.spins must be character strings.');
end
if ~ismember(parameters.spins{1},spin_system.comp.isotopes)
    error('the isotope in parameters.spins is not present in the system.');
end
if size(spectrum,1)~=parameters.zerofill
    error('the number of rows in spectrum must be equal to parameters.zerofill');
end
end

% The only sin on earth is to do things badly.
%
% Ayn Rand

