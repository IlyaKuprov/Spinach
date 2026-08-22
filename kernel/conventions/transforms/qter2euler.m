% Converts a unit quaternion in the active convention into Euler
% angles (ZYZ active convention), matching euler2dcm.m function.
% Syntax:
%
%              [alpha,beta,gamma]=qter2euler(q)
%
% Parameters:
%
%    q - structure with four fields q.u, q.i, q.j, q.k giving
%        the four components of the quaternion; each field may
%        be a column vector, in which case the conversion is
%        performed elementwise
%
% Outputs:
%
%    alpha,beta,gamma - Euler angles in radians (ZYZ active
%                       convention), same shape as the qua-
%                       ternion component fields
%
% Note: Euler angles are not unique; the angles returned sa-
%       tisfy euler2dcm(alpha,beta,gamma)=qter2dcm(q) with
%       beta in the [0,pi] interval.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=qter2euler.m>

function [alpha,beta,gamma]=qter2euler(q)

% Check consistency
grumble(q);

% Normalise the quaternions
qnorm=sqrt(q.u.^2+q.i.^2+q.j.^2+q.k.^2);
q.u=q.u./qnorm; q.i=q.i./qnorm;
q.j=q.j./qnorm; q.k=q.k./qnorm;

% Compute the angle sums and differences
sum_ag=2*atan2(q.k,q.u); dif_ga=2*atan2(q.i,q.j);

% Compute the Euler angles
beta=2*atan2(sqrt(q.i.^2+q.j.^2),sqrt(q.u.^2+q.k.^2));
alpha=(sum_ag-dif_ga)/2; gamma=(sum_ag+dif_ga)/2;

end

% Consistency enforcement
function grumble(q)
if ~all(isfield(q,{'i','j','k','u'}))
    error('quaternion data structure must contain u, i, j, and k fields.');
end
if (~isnumeric(q.u))||(~isreal(q.u))||(~iscolumn(q.u))||...
   (~isnumeric(q.i))||(~isreal(q.i))||(~iscolumn(q.i))||...
   (~isnumeric(q.j))||(~isreal(q.j))||(~iscolumn(q.j))||...
   (~isnumeric(q.k))||(~isreal(q.k))||(~iscolumn(q.k))
    error('quaternion components must be real scalars or column vectors.');
end
if (numel(q.u)~=numel(q.i))||(numel(q.i)~=numel(q.j))||(numel(q.j)~=numel(q.k))
    error('quaternion component fields must have the same number of elements.');
end
if any(sqrt(q.u.^2+q.i.^2+q.j.^2+q.k.^2)<sqrt(eps()))
    error('quaternion norms must be significantly non-zero.');
end
end

