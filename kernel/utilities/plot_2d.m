% Contour plotting utility with non-linear adaptive contour spacing. The
% function is useful for NMR data where small cross-peaks must be adequa-
% tely contoured next to large diagonal peaks. Syntax:
%
%  plot_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs)
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
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=plot_2d.m>

function [axis_f1,axis_f2,spectrum]=plot_2d(spin_system,spectrum,...
                                    parameters,ncont,delta,k,ncol,m,signs)

% Set common defaults
parameters=defaults(spin_system,parameters);

% Check consistency
grumble(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs);

% Inform the user
report(spin_system,'plotting...');

% If a complex spectrum is received, plot both components
if nnz(imag(spectrum))>0
    
    % Recursively plot the real part
    subplot(1,2,1); plot_2d(spin_system,real(spectrum),parameters,ncont,delta,k,ncol,m,signs);
    ktitle('Real part of the complex spectrum.');
    
    % Recursively plot the imaginary part
    subplot(1,2,2); plot_2d(spin_system,imag(spectrum),parameters,ncont,delta,k,ncol,m,signs);
    ktitle('Imaginary part of the complex spectrum.'); return
    
end

% Determine data extents and get contour levels
smax=max(spectrum,[],'all'); smin=min(spectrum,[],'all');
[contours,...
 positive_contours,...
 negative_contours]=contspacing(smax,smin,delta,k,signs,ncont);

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
axis_f1=ft_axis(parameters.offset(1),parameters.sweep(1),size(spectrum,2));
axis_f2=ft_axis(parameters.offset(2),parameters.sweep(2),size(spectrum,1));

% Convert the units 
switch parameters.axis_units
    case 'ppm'
        axis_f1=1000000*(2*pi)*axis_f1/(spin(parameters.spins{1})*spin_system.inter.magnet);
        axis_f2=1000000*(2*pi)*axis_f2/(spin(parameters.spins{2})*spin_system.inter.magnet);
        axis_f1_label=['F1: ' parameters.spins{1} ' chemical shift / ppm'];
        axis_f2_label=['F2: ' parameters.spins{2} ' chemical shift / ppm'];
    case 'Gauss'
        axis_f1=10000*(spin_system.inter.magnet-2*pi*axis_f1/spin('E'));
        axis_f2=10000*(spin_system.inter.magnet-2*pi*axis_f2/spin('E'));
        axis_f1_label='F1: magnetic induction / Gauss';
        axis_f2_label='F2: magnetic induction / Gauss';
    case 'Hz'
        axis_f1=1*axis_f1+0; axis_f2=1*axis_f2+0;
        axis_f1_label=['F1: ' parameters.spins{1} ' linear frequency / Hz'];
        axis_f2_label=['F2: ' parameters.spins{2} ' linear frequency / Hz'];
    case 'kHz'
        axis_f1=0.001*axis_f1+0; axis_f2=0.001*axis_f2+0;
        axis_f1_label=['F1: ' parameters.spins{1} ' linear frequency / kHz'];
        axis_f2_label=['F2: ' parameters.spins{2} ' linear frequency / kHz'];
    case 'MHz'
        axis_f1=0.000001*axis_f1+0; axis_f2=0.000001*axis_f2+0;
        axis_f1_label=['F1: ' parameters.spins{1} ' linear frequency / MHz'];
        axis_f2_label=['F2: ' parameters.spins{2} ' linear frequency / MHz'];
    case 'points'
        axis_f1=1:size(spectrum,2); axis_f2=1:size(spectrum,1);
        axis_f1_label=['F1: ' parameters.spins{1} ' linear frequency / points'];
        axis_f2_label=['F2: ' parameters.spins{2} ' linear frequency / points'];
    otherwise
        error('unknown axis units.');
end

% Plot the spectrum
spectrum=transpose(spectrum);
contour(axis_f2,axis_f1,spectrum,contours);
box on; kgrid; axis square;

% Invert the axes
set(gca,'XDir','reverse','YDir','reverse');

% Label the axes
kxlabel(axis_f2_label); kylabel(axis_f1_label);

% Colour the contours
if any(positive_contours)&&any(negative_contours)
    plot_range=max(abs(positive_contours))+max(abs(negative_contours));
    nredcont=ceil(ncol*max(abs(positive_contours))/plot_range);
    nbluecont=ceil(ncol*max(abs(negative_contours))/plot_range);
elseif any(positive_contours)
    plot_range=max(abs(positive_contours)); nbluecont=0;
    nredcont=ceil(ncol*max(abs(positive_contours))/plot_range);
elseif any(negative_contours)
    plot_range=max(abs(negative_contours)); nredcont=0;
    nbluecont=ceil(ncol*max(abs(negative_contours))/plot_range);
else
    error('spectrum contouring produced no contours.');
end
colors=0.9*(1-[linspace(1,0,nbluecont)' linspace(1,0,nbluecont)' linspace(0,0,nbluecont)';
               linspace(0,0,nredcont)'  linspace(0,1,nredcont)'  linspace(0,1,nredcont)']).^m;
colormap(gca,colors);

% Draw the color bar
if ~ismember('colorbar',spin_system.sys.disable)
    colorbar(gca,'southoutside');
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
function grumble(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs) %#ok<INUSL>
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

% After all, every murderer when he kills runs the risk of the most
% dreadful of deaths, whereas those who kill him risk nothing except
% promotion.
%
% Albert Camus 

