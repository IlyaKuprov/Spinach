% Converts a direction cosine matrix in the active convention of
% euler2dcm.m function into a unit quaternion. Syntax:
%
%                      q=dcm2qter(dcm)
%
% Parameters:
%
%    dcm - directional cosine matrix, a 3x3 orthogonal
%          matrix with unit determinant
%
% Outputs:
%
%    q - structure with four scalar fields q.u, q.i, q.j,
%        q.k giving the four components of the quaternion,
%        normalised to q.u greater than or equal to zero
%
% Note: quaternions double-cover rotations; of the two candi-
%       dates q and -q this function returns the one with the
%       non-negative scalar part.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=dcm2qter.m>

function q=dcm2qter(dcm)

% Check consistency
grumble(dcm);

% Shepperd invariants of the four candidate pivots
piv=[1+dcm(1,1)+dcm(2,2)+dcm(3,3);
     1+dcm(1,1)-dcm(2,2)-dcm(3,3);
     1-dcm(1,1)+dcm(2,2)-dcm(3,3);
     1-dcm(1,1)-dcm(2,2)+dcm(3,3)];

% Convert using the best-conditioned pivot
[~,best]=max(piv);
switch best
    case 1
        q.u=sqrt(piv(1))/2;
        q.i=(dcm(3,2)-dcm(2,3))/(4*q.u);
        q.j=(dcm(1,3)-dcm(3,1))/(4*q.u);
        q.k=(dcm(2,1)-dcm(1,2))/(4*q.u);
    case 2
        q.i=sqrt(piv(2))/2;
        q.u=(dcm(3,2)-dcm(2,3))/(4*q.i);
        q.j=(dcm(1,2)+dcm(2,1))/(4*q.i);
        q.k=(dcm(1,3)+dcm(3,1))/(4*q.i);
    case 3
        q.j=sqrt(piv(3))/2;
        q.u=(dcm(1,3)-dcm(3,1))/(4*q.j);
        q.i=(dcm(1,2)+dcm(2,1))/(4*q.j);
        q.k=(dcm(2,3)+dcm(3,2))/(4*q.j);
    case 4
        q.k=sqrt(piv(4))/2;
        q.u=(dcm(2,1)-dcm(1,2))/(4*q.k);
        q.i=(dcm(1,3)+dcm(3,1))/(4*q.k);
        q.j=(dcm(2,3)+dcm(3,2))/(4*q.k);
end

% Resolve the double cover towards a non-negative scalar part
if q.u<0
    q.u=-q.u; q.i=-q.i; q.j=-q.j; q.k=-q.k;
end

% Normalise the quaternion
qnorm=norm([q.u q.i q.j q.k],2);
q.u=q.u/qnorm; q.i=q.i/qnorm;
q.j=q.j/qnorm; q.k=q.k/qnorm;

end

% Consistency enforcement
function grumble(dcm)
if (~isnumeric(dcm))||(~isreal(dcm))||(~isequal(size(dcm),[3 3]))
    error('dcm must be a real 3x3 matrix.');
end
if norm(dcm'*dcm-eye(3),'fro')>1e-6
    error('dcm must be an orthogonal matrix.');
end
if abs(det(dcm)-1)>1e-6
    error('dcm must have a unit determinant.');
end
end

