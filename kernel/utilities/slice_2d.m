% Contour plotting utility with non-linear adaptive contour spacing and 
% 1D slice extraction using mouse. Syntax:
%
% slice_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)
%
% Parameters:
%
%     spectrum - a real matrix containing the 2D NMR spectrum
%
%     parameters.sweep       -  one or two sweep widths, Hz
%
%     parameters.spins       -  cell array with one ot two character
%                               strings specifying the working spins
%
%     parameters.offset      -  one or two transmitter offsets, Hz
%
%     parameters.axis_units  -  axis units ('ppm','Hz','Gauss')
%
%     ncont     - the number of contours, a reasonable value is 20
%
%     delta     - minimum and maximum elevation (as a fraction of the
%                 total intensity) of the contours above the baseline.
%                 A good starting value is [0.02 0.2 0.02 0.2]. The
%                 first pair of numbers refers to the positive conto-
%                 urs and the second pair to the negative ones.
%
%     k    - a coefficient that controls the curvature of the contour
%            spacing function: k=1 corresponds to linear spacing and
%            k>1 bends the spacing curve to increase the sampling den-
%            sity near the baseline. A reasonable value is 2.
%
%     ncol - number of colours in the colour map; around 256 is fine
%
%     m    - the curvature of the colour map: m=1 corresponds to a li-
%            near colour ramp into the red for positive contours, and
%            into the blue for negative contours. A reasonable value
%            for high-contrast plotting is 6.
%
%     signs   - can be set to 'positive', 'negative' or 'both' - this
%               will cause the corresponding contours to be plotted.
%
% Outputs:
%
%     this function creates a figure
%
% Note: the following functions are used to compute contour levels:
%
%  cont_levs_pos=delta(2)*smax*linspace(0,1,ncont).^k+smax*delta(1);
%  cont_levs_neg=delta(2)*smin*linspace(0,1,ncont).^k+smin*delta(1);
%
% where smin and smax are computed from the spectrum matrix. 
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=slice_2d.m>

function slice_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)

% Set common defaults
parameters=defaults(spin_system,parameters);

% Check consistency
grumble(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs);

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

% Do contour plotting
loc=get(0,'defaultfigureposition'); subplot(1,3,1); 
set(gcf,'Position',[loc(1:2) 3*loc(3) loc(4)]);
[f2,f1,S]=plot_2d(spin_system,spectrum,parameters,...
                  ncont,delta,k,ncol,m,signs);
title('2D spectrum'); drawnow();

% Switch off performance warning
warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId');

% Get the spectrum extents
spec_min=min(spectrum(:));
spec_max=max(spectrum(:));

% Enter the slicing loop
while true()

    % Record pointer position
    subplot(1,3,1); [x,y]=ginput(1);
    
    % Compute axis grids
    [F1,F2]=ndgrid(f1,f2);
    
    % Create the interpolant
    S_int=griddedInterpolant(F1,F2,transpose(S),'spline');
    
    % Compute the traces
    trace_f1=S_int(ones(size(f2))*x,f2);
    trace_f2=S_int(f1,ones(size(f1))*y);
    
    % Call the 1D plotting routine for F1
    parameters_f1=parameters;
    parameters_f1.spins=parameters.spins(1);
    parameters_f1.offset=parameters.offset(1);
    parameters_f1.sweep=parameters.sweep(1);
    parameters_f1.zerofill=parameters.zerofill(1);
    subplot(1,3,2); plot_1d(spin_system,trace_f1,parameters_f1);
    set(gca,'YLim',[spec_min spec_max]); title('F1 slice'); drawnow();
    
    % Call the 1D plotting routine for F2
    parameters_f2=parameters;
    if numel(parameters.spins)==2
        parameters_f2.spins=parameters.spins(2);
    end
    if numel(parameters.offset)==2
        parameters_f2.offset=parameters.offset(2);
    end
    parameters_f2.sweep=parameters.sweep(2);
    parameters_f2.zerofill=parameters.zerofill(2);
    subplot(1,3,3); plot_1d(spin_system,trace_f2,parameters_f2);
    set(gca,'YLim',[spec_min spec_max]); title('F2 slice'); drawnow();
    
end

end

% Default parameters
function parameters=defaults(spin_system,parameters)
if ~isfield(parameters,'offset')
    report(spin_system,'parameters.offset field not set, assuming zero offsets.');
    parameters.offset=zeros(size(parameters.spins));
end
if ~isfield(parameters,'axis_units')
    report(spin_system,'parameters.axis_units field not set, assuming ppm.');
    parameters.axis_units='ppm';
end
end

% Consistency enforcement
function grumble(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)
if (~isnumeric(spectrum))||(~ismatrix(spectrum))
    error('spectrum must be a matrix.');
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
if ~isfield(parameters,'axis_units')
    error('axis units must be specified in parameters.axis_units variable.');
end
if ~ischar(parameters.axis_units)
    error('parameters.axis_units must be a character string.');
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
    if ~ischar(parameters.spins{n})
        error('elements of parameters.spins must be character strings.');
    end
    if ~ismember(parameters.spins{n},spin_system.comp.isotopes)
        error('one of the isotopes in parameters.spins is not present in the system.');
    end
end
if (~isnumeric(ncont))||(~isscalar(ncont))||(~isreal(ncont))||(ncont<1)||(mod(ncont,1)~=0)
    error('ncont parameter must be a positive integer.');
end
if (~isnumeric(delta))||(numel(delta)~=4)||(~isreal(delta))||any(delta>1)||any(delta<0)
    error('delta parameter must be a vector with four real elements between 0 and 1.');
end
if (~isnumeric(k))||(~isscalar(k))||(~isreal(k))||(k<1)||(mod(k,1)~=0)
    error('k parameter must be a positive integer.');
end
if (~isnumeric(ncol))||(~isscalar(ncol))||(~isreal(ncol))||(ncol<1)||(mod(ncol,1)~=0)
    error('ncol parameter must be a positive integer.');
end
if (~isnumeric(m))||(~isscalar(m))||(~isreal(m))||(m<1)||(mod(m,1)~=0)
    error('m parameter must be a positive integer.');
end
if ~ischar(signs)
    error('signs parameter must be a character string.');
end
end

% There's a great deal of difference between an eager man who wants
% to read a book and a tired man who wants a book to read.
%
% Gilbert K. Chesterton

