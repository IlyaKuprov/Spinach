% Principal-axis frame from x and z axes. Syntax:
%
%          frame=diamond_frame_xz(xaxis,zaxis)
%
% This internal helper returns an orthonormal right-handed frame from
% specified x and z axes.
%
% Parameters:
%
%    xaxis  - desired x axis
%
%    zaxis  - desired z axis before orthogonalisation
%
% Outputs:
%
%    frame  - orthonormal right-handed frame
%
% <https://spindynamics.org/wiki/index.php?title=diamond_frame_xz.m>

function frame=diamond_frame_xz(xaxis,zaxis)

% Check consistency
grumble(xaxis,zaxis);

% Normalise the x axis
xaxis=xaxis(:)/norm(xaxis,2);

% Orthogonalise the z axis against x
zaxis=zaxis(:)-xaxis*dot(xaxis,zaxis(:));
zaxis=zaxis/norm(zaxis,2);

% Complete the frame
yaxis=cross(zaxis,xaxis);
frame=diamond_frame_xyz(xaxis,yaxis,zaxis);

end

% Consistency enforcement
function grumble(xaxis,zaxis)
if(~isnumeric(xaxis)||~isreal(xaxis)||numel(xaxis)~=3||norm(xaxis,2)==0)
    error('xaxis must be a non-zero real three-element array.');
end
if(~isnumeric(zaxis)||~isreal(zaxis)||numel(zaxis)~=3||norm(zaxis,2)==0)
    error('zaxis must be a non-zero real three-element array.');
end
end

% Two axes are enough only after orthogonalisation.

