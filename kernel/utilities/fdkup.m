% Returns a finite difference representation of the Kuprov operator:
%
%                 K[rho]=-(1/3)*Trace(Hessian[rho]*chi)
%
% with the number of stencil points in the finite difference approxi-
% mation specified by user. The resulting operator is a sparse matrix
% designed to act on the vectorisation of rho. The dimensions of rho
% are assumed to be ordered as [X Y Z]. For further information, see
% http://dx.doi.org/10.1039/C4CP03106G. Syntax:
%
%                  K=fdkup(npoints,extents,chi,nstenc)
%
% The following parameters are needed:
% 
%     npoints -  a three-element vector specifying the dimensions
%                of the 3D cube of data that the operator will be
%                acting on, in Angstroms. The dimensions are assu-
%                med to be ordered as [X Y Z].
%
%     chi     -  the electron magnetic susceptibility tensor in
%                cubic Angstroms, a symmetric 3x3 matrix.
%
%     extents -  a three-element vector specifying axis extents
%                in Angstroms. The dimensions are assumed to be
%                ordered as [X Y Z].
%
%     nstenc  -  number of finite-difference stencil points for
%                the finite-difference approximation. Periodic 
%                boundary conditions are used.
%
% Outputs:
%
%     K       -  a sparse matrix designed to act on the vectori-
%                zation of the array. The dimensions are assumed
%                to be ordered as [X Y Z].
%
% gareth.charnock@oerc.ox.ac.uk
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=fdkup.m>

function K=fdkup(npoints,extents,chi,nstenc)

% Check consistency
grumble(npoints,extents,chi,nstenc);

% Build second derivative operators
d2_dzdz=kron(kron(fdmat(npoints(3),nstenc,2),speye(npoints(2))),speye(npoints(1)));
d2_dydy=kron(kron(speye(npoints(3)),fdmat(npoints(2),nstenc,2)),speye(npoints(1)));
d2_dxdx=kron(kron(speye(npoints(3)),speye(npoints(2))),fdmat(npoints(1),nstenc,2));
d2_dydx=kron(kron(speye(npoints(3)),fdmat(npoints(2),nstenc,1)),fdmat(npoints(1),nstenc,1));
d2_dzdx=kron(kron(fdmat(npoints(3),nstenc,1),speye(npoints(2))),fdmat(npoints(1),nstenc,1));
d2_dzdy=kron(kron(fdmat(npoints(3),nstenc,1),fdmat(npoints(2),nstenc,1)),speye(npoints(1)));
d2_dxdz=d2_dzdx; d2_dydz=d2_dzdy; d2_dxdy=d2_dydx;

% Normalise second derivative operators
d2_dxdx=(npoints(1)/extents(1))*(npoints(1)/extents(1))*d2_dxdx;
d2_dxdy=(npoints(1)/extents(1))*(npoints(2)/extents(2))*d2_dxdy;
d2_dxdz=(npoints(1)/extents(1))*(npoints(3)/extents(3))*d2_dxdz;
d2_dydx=(npoints(2)/extents(2))*(npoints(1)/extents(1))*d2_dydx;
d2_dydy=(npoints(2)/extents(2))*(npoints(2)/extents(2))*d2_dydy;
d2_dydz=(npoints(2)/extents(2))*(npoints(3)/extents(3))*d2_dydz;
d2_dzdx=(npoints(3)/extents(3))*(npoints(1)/extents(1))*d2_dzdx;
d2_dzdy=(npoints(3)/extents(3))*(npoints(2)/extents(2))*d2_dzdy;
d2_dzdz=(npoints(3)/extents(3))*(npoints(3)/extents(3))*d2_dzdz;

% Form the Kuprov operator
K=-(1/3)*(chi(1,1)*d2_dxdx+chi(1,2)*d2_dxdy+chi(1,3)*d2_dxdz+...
          chi(2,1)*d2_dydx+chi(2,2)*d2_dydy+chi(2,3)*d2_dydz+...
          chi(3,1)*d2_dzdx+chi(3,2)*d2_dzdy+chi(3,3)*d2_dzdz);

end

% Consistency enforcement
function grumble(npoints,extents,chi,nstenc)
if (~isnumeric(npoints))||(numel(npoints)~=3)||(~isreal(npoints))||...
   (any(npoints<1))||any(mod(npoints,1)~=0)
    error('dims must be a three-element vector of positive integers.');
end
if (~isnumeric(extents))||(numel(extents)~=3)||(~isreal(extents))||(any(extents<=0))
    error('extents must be a three-element vector of positive real numbers.');
end
if any(npoints<nstenc)
    error('array dimension is not big enough for the finite difference stencil specified.');
end
if (mod(nstenc,1)~=0)||(mod(nstenc,2)~=1)||(nstenc<3)
    error('the number of stencil points must be an odd integer greater than 3.');
end
if (~isnumeric(chi))||(any(size(chi)~=3))||(~isreal(chi))||(norm(chi-chi',1)>1e-10)
    error('chi must be a real symmetric 3x3 matrix.');
end
end

% The power of accurate observation is commonly called
% cynicism by those who have not got it.
%
% George Bernard Shaw

