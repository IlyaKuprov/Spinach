% Converts a unit quaternion into a direction cosine matrix in
% the active convention, matching the one used by euler2dcm.m
% function. Syntax:
%
%                      dcm=qter2dcm(q)
%
% Parameters:
%
%    q - structure with four scalar fields q.u, q.i, q.j,
%        q.k giving the four components of the quaternion
%
% Outputs:
%
%    dcm - directional cosine matrix
%
% Note: the resulting rotation matrix is to be used as follows:
%
%         v=R*v    (for 3x1 vectors)
%         A=R*A*R' (for 3x3 interaction tensors)
%
% Note: Matlab's Aerospace Toolbox quat2dcm() returns the
%       transpose of this matrix for the same quaternion.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=qter2dcm.m>

function dcm=qter2dcm(q)

% Check consistency
grumble(q);

% Normalise the quaternion
qnorm=norm([q.u q.i q.j q.k],2);
q.u=q.u/qnorm; q.i=q.i/qnorm;
q.j=q.j/qnorm; q.k=q.k/qnorm;

% Compute the DCM
dcm=[1-2*(q.j^2+q.k^2)    2*(q.i*q.j-q.u*q.k)  2*(q.i*q.k+q.u*q.j);
     2*(q.i*q.j+q.u*q.k)  1-2*(q.i^2+q.k^2)    2*(q.j*q.k-q.u*q.i);
     2*(q.i*q.k-q.u*q.j)  2*(q.j*q.k+q.u*q.i)  1-2*(q.i^2+q.j^2)];

end

% Consistency enforcement
function grumble(q)
if ~all(isfield(q,{'i','j','k','u'}))
    error('quaternion data structure must contain u, i, j, and k fields.');
end
if (~isnumeric(q.u))||(~isreal(q.u))||(~isscalar(q.u))||...
   (~isnumeric(q.i))||(~isreal(q.i))||(~isscalar(q.i))||...
   (~isnumeric(q.j))||(~isreal(q.j))||(~isscalar(q.j))||...
   (~isnumeric(q.k))||(~isreal(q.k))||(~isscalar(q.k))
    error('quaternion components must be real scalars.');
end
if norm([q.u q.i q.j q.k],2)<sqrt(eps())
    error('quaternion norm must be significantly non-zero.');
end
end

