% Computes multipole moments from a set of points with user-specified 
% spin populations. Syntax:
%
%                 Ilm=points2mult(xyz,mxyz,rho,L,method)
%
% The multipoles in question are described in Equation 32 of:
%
%                  http://dx.doi.org/10.1039/C6CP05437D
%
% Parameters: 
% 
%    xyz    - coordinates as [x y z] with multiple rows,
%             at which density rho is evaluated, in Angstroms.
% 
%    mxyz   - paramagnetic centre coordinates as [x y z], in 
%             Angstroms.
% 
%    L      - array of ranks of spherical harmonics of the 
%             probability density
% 
%    rho    - column of the densities at the points xyz
% 
%    method - 'points' or 'grid'. If the spin density is supplied 
%             as Mulliken spin populations at individual nuclei,
%             choose 'points'; if the spin density is supplied as
%             a probability on a uniform cubic grid obtained from
%             ndgrid() function and vectorised, use 'grid'.
%
% Output:
% 
%     Ilm  - multipole moments of the probability density 
%
% i.kuprov@soton.ac.uk
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=points2mult.m>

function Ilm=points2mult(xyz,mxyz,rho,L,method)

% Check consistency
grumble(xyz,mxyz,rho,L);

% Put the metal at the origin
xyz=xyz-kron(mxyz,ones(size(xyz,1),1));

% Convert coordinates to spherical
[r,theta,phi]=xyz2sph(xyz(:,1),xyz(:,2),xyz(:,3));

% Compute the coefficients
for n=1:numel(L)
    l=L(n);
    for m=0:l
        if m==0
            Ilm{n}(l+1)=sum(rho.*r.^l.*spher_harmon(l,0,theta,phi));           %#ok<AGROW>
        else
            Ilm{n}( m+l+1)=+real(sum(rho.*r.^l.*spher_harmon(l,m,theta,phi))); %#ok<AGROW>
            Ilm{n}(-m+l+1)=-imag(sum(rho.*r.^l.*spher_harmon(l,m,theta,phi))); %#ok<AGROW>
        end
    end
end

% Process the options
switch method
    
    case 'grid'
        
        % Compute the volume element
        x=unique(xyz(:,1)); dx=(max(x)-min(x))/(numel(x)-1);
        y=unique(xyz(:,2)); dy=(max(y)-min(y))/(numel(y)-1);
        z=unique(xyz(:,3)); dz=(max(z)-min(z))/(numel(z)-1);
        
        % Apply the volume element
        Ilm=Ilm*dx*dy*dz;
        
    case 'points'
        
        % No action necessary
        
end

end

% Consistency enforcement
function grumble(xyz,mxyz,rho,L)
if (~isnumeric(xyz))||(~isreal(xyz))||(size(xyz,2)~=3)
    error('xyz must be an Nx3 array of atomic coordinates.');
end
if (~isnumeric(rho))||(size(rho,2)~=1)||(~isreal(rho))
    error('rho parameter should be a real column vector.');
end
if size(xyz,1)~=size(rho,1)
    error('the number of rows in xyz and rho arguments must be the same.');
end
if (~isnumeric(mxyz))||(size(mxyz,2)~=3)||(size(mxyz,1)~=1)||(~isreal(mxyz))
    error('mxyz parameter should be a real row vector with three elements.');
end
if (~isnumeric(L))||(~isreal(L))
    error('L must be a real vector.');
end
end

% You can fall out with the master; that doesn't matter at all. The 
% three people in a college you must never fall out with are the head
% gardener, the head porter and the head chef.
%
% Rory Sutherland

