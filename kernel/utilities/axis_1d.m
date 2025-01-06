% Generates axis ticks for plotting 1D spectra. Syntax:
%
%             [ax,ax_label]=axis_1d(spin_system,parameters)
%
% Inputs:
%
%    parameters.sweep  -  either a one-element array giving the sweep width
%                         in Hz, or a two-element array giving the spectral
%                         extents in Hz around the offset.
%
%    parameters.zerofill - the number of points in the NMR spectrum after 
%                          zerofilling and Fourier transform
%
%    parameters.offset   - offset of the spectrum centre point relative to
%                          the magnet frequency, Hz
%
% parameters.axis_units  - a character string with the units in which the
%                          axis ticks should be returned: 'ppm', 'Gauss',
%                          'mT', 'Hz', 'kHz', 'MHz','MHz-labframe', 'GHz',
%                          'GHz-labframe', 'gtensor', 'points'
%
%    parameters.spins      - the spin involved, e.g. {'1H'}
%
% Outputs:
%
%    axis                  - a row vector of axis tick values
%
%    ax_label              - axis label for displaying on the plot
%
% Note: magnetic field units use the free electron g-tensor for conversion.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=axis_1d.m>

function [ax,ax_label]=axis_1d(spin_system,parameters)

% Check consistency
grumble(spin_system,parameters);

% Build the axis
if isscalar(parameters.sweep)
    
    % Offset + sweep specification
    ax=ft_axis(parameters.offset,parameters.sweep,parameters.zerofill);
             
    % Report to the user
    report(spin_system,'offset + sweep frequency axis specification:');
    report(spin_system,[num2str(parameters.zerofill) '-point axis from '...
                        num2str(ax(1)) ' to ' num2str(ax(end)) ' Hz']);
    
else
    
    % Two-element sweep specification
    ax=linspace(parameters.sweep(1),parameters.sweep(2),parameters.zerofill);
            
    % Report to the user
    report(spin_system,'two-element sweep frequency axis specification:');
    report(spin_system,[num2str(parameters.zerofill) '-point axis from '...
                        num2str(ax(1)) ' to ' num2str(ax(end)) ' Hz']);
    
end

% Carrier frequency of the spin in Hz
basefrq=-spin(parameters.spins{1})*spin_system.inter.magnet/(2*pi);

% Convert the units if necessary
switch parameters.axis_units
    
    case 'ppm'
        
        % NMR-style ppm scale with respect to the carrier frequency
        ax_label=[parameters.spins{1} ' chemical shift / ppm'];
        ax=-1e6*ax/basefrq;
        
    case 'Gauss'
        
        % Pulsed FFT EPR style Gauss axis for rotating frame simulations
        ax_label='Magnetic induction / Gauss';
        ax=1e4*(spin_system.inter.magnet-2*pi*ax/spin('E'));
        
    case 'mT'
        
        % Pulsed FFT EPR style mT axis for rotating frame simulations
        ax_label='Magnetic induction / mT';
        ax=1e3*(spin_system.inter.magnet-2*pi*ax/spin('E'));
        
    case 'Hz'
        
        % Raw Hz axis for rotating frame simulations
        ax_label=[parameters.spins{1} ' offset frequency / Hz'];
        
    case 'kHz'
        
        % Raw kHz axis for rotating frame simulations
        ax_label=[parameters.spins{1} ' offset frequency / kHz'];
        ax=1e-3*ax;
        
    case 'MHz'
        
        % Raw MHz axis for rotating frame simulations
        ax_label=[parameters.spins{1} ' offset frequency / MHz'];
        ax=1e-6*ax;
        
    case 'MHz-labframe'
        
        % Raw MHz axis for lab frame simulations
        ax_label=[parameters.spins{1} ' frequency / MHz'];
        ax=-1e-6*(ax-basefrq);
        
    case 'GHz'
        
        % Raw GHz axis for rotating frame simulations
        ax_label=[parameters.spins{1} ' offset frequency / GHz'];
        ax=1e-9*ax;
        
    case 'GHz-labframe'
        
        % Raw GHz axis for lab frame simulations
        ax_label=[parameters.spins{1} ' frequency / GHz'];
        ax=-1e-9*(ax-basefrq);
        
    case 'gtensor'
        
        % Raw g-tensor axis for EPR simulations
        ax_label=[parameters.spins{1} ' g-tensor / \mu_B'];
        ax=-spin_system.tols.freeg*(ax-basefrq)/basefrq;
        
    case 'points'
        
        % Raw point count
        ax_label='digitisation points';
        ax=1:parameters.zerofill;
        
    otherwise
        
        % Complain and bomb out
        error('unknown axis units.');
        
end

end

% Consistency enforcement
function grumble(spin_system,parameters)
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
if (~isfield(spin_system,'inter'))||(~isfield(spin_system.inter,'magnet'))||...
   (spin_system.inter.magnet==0)
    if strcmp(parameters.axis_units,'ppm')
        error('ppm scale is not available when magnet field is zero.');
    end
end
end

% The God of the Old Testament is arguably the most unpleasant character in 
% all fiction: jealous and proud of it; a petty, unjust, unforgiving control
% freak; a vindictive, bloodthirsty ethnic cleanser; a misogynistic, homopho-
% bic, racist, infanticidal, genocidal, filicidal, pestilential, megalomania-
% cal, sadomasochistic, capriciously malevolent bully.
%
% Richard Dawkins, "The God Delusion"

