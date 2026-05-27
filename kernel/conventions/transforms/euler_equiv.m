% Checks whether two ZYZ active Euler angle sets specify the same
% rotation. Syntax:
%
%              answer=euler_equiv(eulers_a,eulers_b,tol)
%
% Parameters:
%
%    eulers_a   - first Euler angle set [alpha beta gamma],
%                 radians, ZYZ active convention
%
%    eulers_b   - second Euler angle set [alpha beta gamma],
%                 radians, ZYZ active convention
%
%    tol        - non-negative angular tolerance, radians
%
% Outputs:
%
%    answer     - true if the relative rotation angle between
%                 the two rotations is not greater than tol
%
% Note: Euler angles are not unique, and so this function compares
%       the rotations produced by euler2dcm(), not the angles
%       themselves.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=euler_equiv.m>

function answer=euler_equiv(eulers_a,eulers_b,tol)

% Check consistency
grumble(eulers_a,eulers_b,tol);

% Cast Euler angles to row vectors
eulers_a=eulers_a(:)';
eulers_b=eulers_b(:)';

% Build direction cosine matrices
dcm_a=euler2dcm(eulers_a);
dcm_b=euler2dcm(eulers_b);

% Compute the relative rotation
dcm_rel=dcm_b*dcm_a';

% Compute the sine and cosine of the relative rotation angle
sin_ang=norm([dcm_rel(3,2)-dcm_rel(2,3);...
              dcm_rel(1,3)-dcm_rel(3,1);...
              dcm_rel(2,1)-dcm_rel(1,2)],2)/2;
cos_ang=(trace(dcm_rel)-1)/2;

% Clamp round-off in the cosine
cos_ang=max(min(cos_ang,1),-1);

% Compare the geodesic distance on SO(3)
answer=(atan2(sin_ang,cos_ang)<=tol);

end

% Consistency enforcement
function grumble(eulers_a,eulers_b,tol)
if (~isnumeric(eulers_a))||(~isreal(eulers_a))||(~isvector(eulers_a))||...
   (numel(eulers_a)~=3)||any(~isfinite(eulers_a(:)))
    error('eulers_a must be a real finite vector with three elements.');
end
if (~isnumeric(eulers_b))||(~isreal(eulers_b))||(~isvector(eulers_b))||...
   (numel(eulers_b)~=3)||any(~isfinite(eulers_b(:)))
    error('eulers_b must be a real finite vector with three elements.');
end
if (~isnumeric(tol))||(~isreal(tol))||(numel(tol)~=1)||...
   (~isfinite(tol))||(tol<0)
    error('tol must be a non-negative real finite scalar.');
end
end


