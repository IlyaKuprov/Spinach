% Axis-angle rotation matrix. Syntax:
%
%          R=diamond_rot_axis(axis,angle)
%
% This internal helper returns a rotation matrix for an axis-angle
% rotation in degrees.
%
% Parameters:
%
%    axis   - rotation axis
%
%    angle  - rotation angle, degrees
%
% Outputs:
%
%    R      - rotation matrix
%
% <https://spindynamics.org/wiki/index.php?title=diamond_rot_axis.m>

function R=diamond_rot_axis(axis,angle)

% Check consistency
grumble(axis,angle);

% Normalise the rotation axis
axis=axis(:)/norm(axis,2);

% Build the cross-product matrix
K=[0 -axis(3) axis(2);axis(3) 0 -axis(1);-axis(2) axis(1) 0];

% Apply Rodrigues' rotation formula
R=eye(3)+sind(angle)*K+(1-cosd(angle))*(K*K);

end

% Consistency enforcement
function grumble(axis,angle)
if(~isnumeric(axis)||~isreal(axis)||numel(axis)~=3||norm(axis,2)==0)
    error('axis must be a non-zero real three-element array.');
end
if(~isnumeric(angle)||~isreal(angle)||~isscalar(angle))
    error('angle must be a real scalar.');
end
end

% Small angular corrections should remain explicit rotations.

