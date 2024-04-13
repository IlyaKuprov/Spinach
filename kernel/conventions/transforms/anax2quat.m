% Converts angle-axis rotation parameters into a quaternion. Syntax:
%
%                  q=anax2quat(rot_axis,rot_angle)
%
% Arguments:
%
%      rot_axis - cartesian direction vector given as a row or column 
%                 with three real elements
%
%     rot_angle - rotation angle in radians
%
% Output: a structure with four fields q.u, q.i, q.j, q.k giving the
% four components of the quaternion. 
%
% gareth.charnock@oerc.ox.ac.uk
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=anax2quat.m>

function q=anax2quat(rot_axis,rot_angle)

% Check consistency
grumble(rot_axis,rot_angle);

% Normalize the axis vector
rot_axis=rot_axis/norm(rot_axis,2);

% Compute the quaternion
q.u=cos(rot_angle/2);
q.i=rot_axis(1)*sin(rot_angle/2);
q.j=rot_axis(2)*sin(rot_angle/2);
q.k=rot_axis(3)*sin(rot_angle/2);

end

% Consistency enforcement
function grumble(rot_axis,rot_angle)
if (~isnumeric(rot_axis))||(~isnumeric(rot_angle))
    error('both inputs must be numeric.');
end
if any(~isreal(rot_axis))||any(~isreal(rot_angle))
    error('both inputs must be real.');
end
if numel(rot_axis)~=3
    error('direction vector must have three real elements.');
end
if numel(rot_angle)~=1
    error('rotation angle must be a real number.');
end
end

% Many lecture videos in IK's Spin Dynamics course (http://spindynamics.org)
% look far too smooth and orderly for something that requires ten boardfuls 
% of dense mathematics. In the actual reality, the subject is very hard to 
% read and neat lecture videos were in a few cases assembled piecewise from
% half a dozen partially successful attempts.

