% Converts Euler angles (ZYZ active convention) into a unit
% quaternion in the active convention, matching euler2dcm.m
% function. Syntax:
%
%                q=euler2qter(alpha,beta,gamma)
%
%                            OR
%
%                q=euler2qter([alpha beta gamma])
%
% Parameters:
%
%    alpha,beta,gamma - Euler angles in radians (ZYZ active
%                       convention), scalars or column vec-
%                       tors of equal length
%
% Outputs:
%
%    q - structure with four fields q.u, q.i, q.j, q.k giving
%        the four components of the quaternion; for column
%        vector inputs each field is a column vector
%
% Note: the quaternion returned represents the same rotation
%       as euler2dcm(alpha,beta,gamma); it is converted into
%       that matrix by qter2dcm.m function.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=euler2qter.m>

function q=euler2qter(arg1,arg2,arg3)

% Adapt to the input style
if nargin==1

    % Assume that a single input is a 3-vector
    alp=arg1(1); bet=arg1(2); gam=arg1(3);

elseif nargin==3

    % Assume that three inputs are Euler angles
    alp=arg1; bet=arg2; gam=arg3;

else

    % Bomb out in all other cases
    error('incorrect number of input arguments.');

end

% Check consistency
grumble(alp,bet,gam);

% Compute the quaternion
q.u=cos(bet/2).*cos((alp+gam)/2);
q.i=sin(bet/2).*sin((gam-alp)/2);
q.j=sin(bet/2).*cos((gam-alp)/2);
q.k=cos(bet/2).*sin((alp+gam)/2);

end

% Consistency enforcement
function grumble(alp,bet,gam)
if (~isnumeric(alp))||(~isreal(alp))||(~iscolumn(alp))||...
   (~isnumeric(bet))||(~isreal(bet))||(~iscolumn(bet))||...
   (~isnumeric(gam))||(~isreal(gam))||(~iscolumn(gam))
    error('Euler angles must be real scalars or column vectors.');
end
if (numel(alp)~=numel(bet))||(numel(bet)~=numel(gam))
    error('alpha, beta, and gamma must have the same number of elements.');
end
end

