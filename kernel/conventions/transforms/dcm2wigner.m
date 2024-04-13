% Converts a directional cosine matrix into second-rank Wigner function 
% matrix. Syntax:
%
%                             D=dcm2wigner(dcm)
%
% Parameters:
%
%     dcm     - directional cosine matrix
%
% Outputs:
%
%     D       - matrix of second rank Wigner D functions. Rows
%               and columns are sorted by descending ranks:
%    
%                         [D( 2,2)  ...  D( 2,-2)
%                            ...    ...    ...  
%                          D(-2,2)  ...  D(-2,-2)]
%
% Notes: the resulting Wigner matrix is to be used as v=W*v, where v is
%        a column vector of irreducible spherical tensor coefficients in
%        the following order: T(2,2), T(2,1), T(2,0), T(2,-1), T(2,-2).
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=dcm2wigner.m>

function D=dcm2wigner(dcm)

% Check consistency
grumble(dcm);

% Compute A and B coefficients
A=sqrt(0.5*(+dcm(1,1)+1i*dcm(1,2)-1i*dcm(2,1)+dcm(2,2)));
B=sqrt(0.5*(-dcm(1,1)+1i*dcm(1,2)+1i*dcm(2,1)+dcm(2,2)));

% Verify amplitudes
if abs(A*A'-0.5*(1+dcm(3,3)))+abs(B*B'-0.5*(1-dcm(3,3)))>1e-6
    error('DCM does not pass self-consistency check on amplitudes.');
end

% Verify phases
if abs(A*B+0.5*(dcm(1,3)-1i*dcm(2,3)))+abs(A*B'-0.5*(dcm(3,1)+1i*dcm(3,2)))>1e-6
    A=-A;
end
if abs(A*B+0.5*(dcm(1,3)-1i*dcm(2,3)))+abs(A*B'-0.5*(dcm(3,1)+1i*dcm(3,2)))>1e-6
    error('DCM does not pass self-consistency check on phases.');
end

% Compute Wigner matrix
Z=A*A'-B*B';
D=[ A^4                 2*A^3*B           sqrt(6)*A^2*B^2      2*A*B^3           B^4             
   -2*A^3*B'            A^2*(2*Z-1)       sqrt(6)*A*B*Z        B^2*(2*Z+1)       2*A'*B^3         
    sqrt(6)*A^2*B'^2   -sqrt(6)*A*B'*Z    0.5*(3*Z^2-1)        sqrt(6)*A'*B*Z    sqrt(6)*A'^2*B^2 
   -2*A*B'^3            B'^2*(2*Z+1)     -sqrt(6)*A'*B'*Z      A'^2*(2*Z-1)      2*A'^3*B         
    B'^4               -2*A'*B'^3         sqrt(6)*A'^2*B'^2   -2*A'^3*B'         A'^4             ];

end

% Consistency enforcement
function grumble(dcm)
if (~isnumeric(dcm))||(~isreal(dcm))||(~all(size(dcm)==[3 3]))
    error('DCM must be a real 3x3 matrix.');
end
if norm(dcm'*dcm-eye(3),1)>1e-6
    warning('DCM is not orthogonal to 1e-6 tolerance - conversion accuracy not guaranteed.');
end
if norm(dcm'*dcm-eye(3),1)>1e-2
    error('DCM is not orthogonal to 1e-2 tolerance, cannot proceed with conversion.');
end
if abs(det(dcm)-1)>1e-6
    warning('DCM determinant is not unit to 1e-6 tolerance - conversion accuracy not guaranteed.');
end
if abs(det(dcm)-1)>1e-2
    error('DCM determinant is not unit to 1e-2 tolerance, cannot proceed with conversion.');
end
end

% The hallmark of a second rater is resentment for 
% another man's achievement.
%
% Ayn Rand, "Atlas Shrugged"

