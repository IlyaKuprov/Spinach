% Principal-axis frame for Nadolinny cobalt tables. Syntax:
%
%          frame=diamond_frame_cobalt(alpha)
%
% This internal helper returns the frame convention used for cobalt
% alpha-angle tables.
%
% Parameters:
%
%    alpha  - table angle in degrees
%
% Outputs:
%
%    frame  - principal-axis frame
%
% <https://spindynamics.org/wiki/index.php?title=diamond_frame_cobalt.m>

function frame=diamond_frame_cobalt(alpha)

% Check consistency
grumble(alpha);

% Set the fixed table axes
yaxis=[0;1;1]/sqrt(2);
xaxis=[1;0;0];
zbase=cross(xaxis,yaxis);

% Rotate the x axis by the table angle
xaxis=cosd(alpha)*xaxis+sind(alpha)*zbase;

% Complete the frame
zaxis=cross(xaxis,yaxis);
frame=diamond_frame_xyz(xaxis,yaxis,zaxis);

end

% Consistency enforcement
function grumble(alpha)
if(~isnumeric(alpha)||~isreal(alpha)||~isscalar(alpha))
    error('alpha must be a real scalar.');
end
end

% Cobalt tables have their own compass.

