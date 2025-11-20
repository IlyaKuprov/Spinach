% Bloch sphere plot of a 3D trajectory with linear inter-
% polation between points. Syntax:
%
%              bloch_sph_plot(x,y,z,varargin)
%
% Parameters:
%
%   x,y,z    - arrays of equal size containing
%              trajectory or trajectories (rows)
%
%   varargin - further arguments to be supplied 
%              to the plot3() functon of Matlab
%
% Output:
%
%   this function produces a figure
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=bloch_sph_plot.m>

function bloch_sph_plot(x,y,z,varargin)

% Kill stupid ass figure defaults in R2025a and later 
set(groot,'defaultFigurePosition',[680 458 560 420]); 
set(groot,'defaultFigureWindowStyle','normal'); 
set(groot,'defaultFigureMenuBar','figure'); 
set(groot,'defaultFigureToolbar','figure'); 

% Reference unit sphere
figure(); [X,Y,Z]=sphere; 
surf(X,Y,Z,'FaceAlpha',0.25,'EdgeAlpha',0.25);
colormap bone; box on; axis square; hold on;
set(gca,'Projection','perspective'); kgrid;

% Trajectories
plot3(x,y,z,varargin{:});

end

% The Proton VPN servers for Afghanistan running at 
% over 80% load was not something that we expected
% to see when we set them up last year. 
%
% David Peterson, after the UK started
% requiring proof of age on porn sites
% in August 2025

