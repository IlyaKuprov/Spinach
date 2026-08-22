% Converts angle-axis rotation parameters to a direction
% cosine matrix in the active convention, matching the one
% used by euler2dcm.m function. Angle should be in radians,
% axis is normalized by the function. Syntax:
%
%              dcm=anax2dcm(rot_axis,rot_angle)
%
% Parameters:
%
%      rot_axis - cartesian direction vector given as
%                 a row or column with three real ele-
%                 ments
%
%     rot_angle - rotation angle in radians
%
% Outputs:
%
%           dcm - directional cosine matrix
%
% Note: the resulting rotation matrix is to be used as follows:
%
%         v=R*v    (for 3x1 vectors)
%         A=R*A*R' (for 3x3 interaction tensors)
%
% Note: Matlab's Aerospace Toolbox quat2dcm() returns the
%       transpose of this matrix for the same rotation.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=anax2dcm.m>

function dcm=anax2dcm(rot_axis,rot_angle)

% Check consistency
grumble(rot_axis,rot_angle);

% Normalize the axis
rot_axis=rot_axis(:)/norm(rot_axis(:),2);

% Compute the DCM
dcm=eye(3)+sin(rot_angle)*[ 0           -rot_axis(3)  rot_axis(2);
                            rot_axis(3)  0           -rot_axis(1);
                           -rot_axis(2)  rot_axis(1)  0          ]+...
                    (1-cos(rot_angle))*(rot_axis*rot_axis'-eye(3));

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
if norm(rot_axis,2)==0
    error('direction vector must be non-zero.');
end
if numel(rot_angle)~=1
    error('rotation angle must be a real number.');
end
end

% "The brightest flame casts the darkest shadow."
%
% George R.R. Martin

