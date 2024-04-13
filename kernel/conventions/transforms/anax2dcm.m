% Converts angle-axis rotation parameters to directional
% cosine matrix. Angle should be in radians, axis is nor-
% malized by the function. Syntax:
%
%              dcm=anax2dcm(rot_axis,rot_angle)
%
% Arguments:
%
%      rot_axis - cartesian direction vector given as 
%                 a row or column with three real ele-
%                 ments
%
%     rot_angle - rotation angle in radians
%
% Ouputs:
%
%           dcm - directional cosine matrix
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=anax2dcm.m>

function dcm=anax2dcm(rot_axis,rot_angle)

% Check consistency
grumble(rot_axis,rot_angle);

% Normalize the axis
rot_axis=rot_axis(:)/norm(rot_axis(:),2);

% Compute the DCM
dcm=eye(3)-sin(rot_angle)*[ 0           -rot_axis(3)  rot_axis(2);
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
if numel(rot_angle)~=1
    error('rotation angle must be a real number.');
end
end

% "The brightest flame casts the darkest shadow."
%
% George R.R. Martin

