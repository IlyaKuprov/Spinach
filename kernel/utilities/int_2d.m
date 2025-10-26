% Contour plotting utility with non-linear adaptive contour spacing and 
% 2D integration using mouse or an interval file. Syntax:
%
%              int_2d(spin_system,spectrum,parameters,ncont,...
%                     delta,k,ncol,m,signs,filename)
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
%     filename  - the name of the interval file. If this file does not
%                 exist, the integration intervals are saved into this
%                 file. If it does exist, the function prints the in-
%                 tegrals over the ontervals found in this file.
%
% Outputs:
%
%     this function creates a figure and writes the interval file
%
% Note: the following functions are used to compute contour levels:
%
%  cont_levs_pos=delta(2)*smax*linspace(0,1,ncont).^k+smax*delta(1);
%  cont_levs_neg=delta(2)*smin*linspace(0,1,ncont).^k+smin*delta(1);
%
% where smin and smax are computed from the spectrum matrix. 
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=int_2d.m>

function int_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs,filename)

% Check consistency
grumble(filename);

% Do contour plotting (a detailed grumbler inside)
[f2,f1,S]=plot_2d(spin_system,spectrum,parameters,ncont,delta,k,ncol,m,signs);

% Switch off performance warning
warning('off','MATLAB:griddedInterpolant:MeshgridEval2DWarnId');

% Check if the file exists
if ~exist(filename,'file')
    
    % Get the counter going
    n=1;
    
    % Report to the user
    disp('Interactive integration, define the box by clicking on its opposite corners.');
    disp(['Saving integration ranges into ' filename ' ...']);
    
    % Proceed with interactive integration
    while true()
        
        % Get mouse input
        [ranges{n,1},ranges{n,2}]=ginput(2); %#ok<AGROW>
        
        % Write the file
        save(filename,'ranges'); drawnow;
        
        % Compute axis grids
        [F1,F2]=ndgrid(f1,f2);
        
        % Create the interpolant
        S_int=griddedInterpolant(F1,F2,transpose(S),'spline');
        
        % Make function handle
        fhandle=@(x,y)S_int(x,y);
        
        % Compute the integral
        I=integral2(fhandle,ranges{n,1}(1),ranges{n,1}(2),...
                            ranges{n,2}(1),ranges{n,2}(2),'RelTol',1e-3,'AbsTol',1e-3);
        
        % Report to the user
        disp(['Recorded X range: ' num2str(min(ranges{n,1})) ' to ' num2str(max(ranges{n,1}))]);
        disp(['Recorded Y range: ' num2str(min(ranges{n,2})) ' to ' num2str(max(ranges{n,2}))]);
        disp(['Integral: ' num2str(abs(I))]);
        
        % Increment counter
        n=n+1;
        
    end
    
else
    
    % Load ranges from file
    disp('Found the ranges file, integrals:'); 
    load(filename,'ranges');
    
    % Perform automatic integration
    for n=1:size(ranges,1)
        
        % Compute axis grids
        [F1,F2]=ndgrid(f1,f2);
        
        % Create the interpolant
        S_int=griddedInterpolant(F1,F2,transpose(S),'spline');
        
        % Make function handle
        fhandle=@(x,y)S_int(x,y);
        
        % Compute the integral
        I=integral2(fhandle,ranges{n,1}(1),ranges{n,1}(2),...
                            ranges{n,2}(1),ranges{n,2}(2),'RelTol',1e-3,'AbsTol',1e-3);
        
        % Report to the user
        disp(I);
        
    end
    
end

end

% Consistency enforcement
function grumble(filename)
if ~ischar(filename)
    error('filename argument must be a character string.');
end
end

% Beware of the person of one book.
%
% Thomas Aquinas

