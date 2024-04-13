% Converts a quaternion representation of a rotation into angle-axis
% rotation parameters. Syntax:
%
%                 [rot_axis,rot_angle]=quat2anax(q)
%
% Parameters:
%
%     q   -   quaternion, a structure with four fields
%             q.u, q.i, q.j, q.k giving the four compo-
%             nents of the quaternion
%
% Outputs:
%
%     rot_axis  - cartesian direction vector as a row 
%                 with three real elements
%
%     rot_angle - rotation angle in radians
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=quat2anax.m>

function [rot_axis,rot_angle]=quat2anax(q)

% Check consistency
grumble(q);

% Normalize the quaternion
qnorm=norm([q.u q.i q.j q.k],2);
q.u=q.u/qnorm; q.i=q.i/qnorm;
q.j=q.j/qnorm; q.k=q.k/qnorm;

% Compute the angle
rot_angle=2*atan2(norm([q.i q.j q.k],2),q.u);

% Compute the axis
if rot_angle==0
    rot_axis=[0 0 1];
else
    rot_axis=[q.i q.j q.k]/norm([q.i q.j q.k],2);
end

end

% Consistency enforcement
function grumble(q)
if ~all(isfield(q,{'i','j','k','u'}))
    error('quaternion data structure must contain u, i, j, and k fields.');
end
if ~isreal([q.u q.i q.j q.k])
    error('quaternion elements must be real.');
end
end

% "Mediocrity" does not mean an average intelligence; it means an
% average intelligence that resents and envies its betters.
%
% Ayn Rand

