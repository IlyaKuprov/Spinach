% Principal-axis frame for Nadolinny alpha tables. Syntax:
%
%          frame=diamond_frame_alpha(alpha)
%
% This internal helper returns the frame convention used for the
% nickel and titanium alpha-angle tables.
%
% Parameters:
%
%    alpha  - table angle in degrees
%
% Outputs:
%
%    frame  - principal-axis frame
%
% <https://spindynamics.org/wiki/index.php?title=diamond_frame_alpha.m>

function frame=diamond_frame_alpha(alpha)

% Check consistency
grumble(alpha);

% Set the fixed table axes
xaxis=[1;-1;0]/sqrt(2);
ybase=[1;1;0]/sqrt(2);
zbase=[0;0;1];

% Rotate the y axis by the table angle
yaxis=cosd(alpha)*ybase+sind(alpha)*zbase;

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

% Tabulated angles should not hide their axis convention.

