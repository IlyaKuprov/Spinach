% Principal-axis frame from a z axis. Syntax:
%
%          frame=diamond_frame_z(zaxis)
%
% This internal helper returns a right-handed principal-axis frame with
% the specified z axis.
%
% Parameters:
%
%    zaxis  - desired z axis
%
% Outputs:
%
%    frame  - orthonormal right-handed frame
%
% <https://spindynamics.org/wiki/index.php?title=diamond_frame_z.m>

function frame=diamond_frame_z(zaxis)

% Check consistency
grumble(zaxis);

% Normalise the specified z axis
zaxis=zaxis(:)/norm(zaxis,2);

% Pick a non-collinear x-axis seed
if abs(dot(zaxis,[0;0;1]))<0.9
    xaxis=cross([0;0;1],zaxis);
else
    xaxis=cross([0;1;0],zaxis);
end

% Complete the right-handed frame
xaxis=xaxis/norm(xaxis,2);
yaxis=cross(zaxis,xaxis);
frame=[xaxis yaxis zaxis];

end

% Consistency enforcement
function grumble(zaxis)
if(~isnumeric(zaxis)||~isreal(zaxis)||numel(zaxis)~=3||norm(zaxis,2)==0)
    error('zaxis must be a non-zero real three-element array.');
end
end

% A symmetry axis still needs a complete frame.

