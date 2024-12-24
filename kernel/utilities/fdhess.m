% Returns the finite-difference Hessian of a 3D array using a finite
% difference scheme with a user-specified number of stencil points and
% a unit grid spacing. The dimensions of the 3D array are assumed to
% be ordered as [X Y Z]. Syntax:
%
%                         H=fdhess(A,nstenc)
%
% Parameters:
%
%     A      - a 3D array with dimensions ordered as [X Y Z]
%
%     nstenc - number of points inthe finite difference stencil,
%              periodic boundary conditions are used
%
% Outputs:
%
%     H      - a 3x3 cell array of 3D arrays ordered in the
%              following way:
%
%                     {d2A_dxdx  d2A_dxdy  d2A_dxdz
%                      d2A_dydx  d2A_dydy  d2A_dydz
%                      d2A_dzdx  d2A_dzdy  d2A_dzdz}
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=fdhess.m>

function H=fdhess(A,nstenc)

% Check consistency
grumble(A,nstenc);

% Compute derivatives
d2A_dzdz=reshape(kron(kron(fdmat(size(A,3),nstenc,2),speye(size(A,2))),speye(size(A,1)))*A(:),size(A));
d2A_dzdy=reshape(kron(kron(fdmat(size(A,3),nstenc,1),fdmat(size(A,2),nstenc,1)),speye(size(A,1)))*A(:),size(A));
d2A_dzdx=reshape(kron(kron(fdmat(size(A,3),nstenc,1),speye(size(A,2))),fdmat(size(A,1),nstenc,1))*A(:),size(A));
d2A_dydz=reshape(kron(kron(fdmat(size(A,3),nstenc,1),fdmat(size(A,2),nstenc,1)),speye(size(A,1)))*A(:),size(A));
d2A_dydy=reshape(kron(kron(speye(size(A,3)),fdmat(size(A,2),nstenc,2)),speye(size(A,1)))*A(:),size(A));
d2A_dydx=reshape(kron(kron(speye(size(A,3)),fdmat(size(A,2),nstenc,1)),fdmat(size(A,1),nstenc,1))*A(:),size(A));
d2A_dxdz=reshape(kron(kron(fdmat(size(A,3),nstenc,1),speye(size(A,2))),fdmat(size(A,1),nstenc,1))*A(:),size(A));
d2A_dxdy=reshape(kron(kron(speye(size(A,3)),fdmat(size(A,2),nstenc,1)),fdmat(size(A,1),nstenc,1))*A(:),size(A));
d2A_dxdx=reshape(kron(kron(speye(size(A,3)),speye(size(A,2))),fdmat(size(A,1),nstenc,2))*A(:),size(A));

% Form the Hessian array
H={d2A_dxdx d2A_dxdy d2A_dxdz; 
   d2A_dydx d2A_dydy d2A_dydz; 
   d2A_dzdx d2A_dzdy d2A_dzdz};

end

% Consistency enforcement
function grumble(A,npoints)
if (~isnumeric(A))||(ndims(A)~=0)
    error('A must be a three-dimensional numeric array.');
end
if any(size(A)<npoints)
    error('the dimension of A is not big enough for the finite difference stencil specified.');
end
if (mod(npoints,1)~=0)||(mod(npoints,2)~=1)||(npoints<3)
    error('the number of stencil points must be an odd integer greater than 3.');
end
end

% Take therefore the talent from him, and give it unto him which
% hath ten talents. For unto every one that hath shall be given,
% and he shall have abundance: but from him that hath not shall
% be taken away even that which he hath. And cast ye the unprofi-
% table servant into outer darkness: there shall be weeping and
% gnashing of teeth.
%
% Matthew 25:28-30

