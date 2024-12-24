% Computes the three-dimensional distribution of pseudocontact shift field
% by solving Kuprov equation for PCS. Syntax:
%
%            [pcs_cube,pcs_vals]=kpcs(probden,chi,ranges,nxyz)
%
% Parameters: 
%
%     probden - electron probability density cube.
%
%     chi     - electron magnetic susceptibility tensor in cubic Angstroms.
%
%     ranges  - Cartesian axis extents for the electron spin density cube
%               as [xmin xmax ymin ymax zmin zmax] in angstroms.
%
%     nxyz    - nuclear coordinates as [x y z] with multiple rows) at which
%               PCS is to be evaluated, in Angstroms.
%
% Output:
% 
%     pcs_vals - pseudocontact shift in ppm at each nucleus.
%
%     pcs_cube - pseudocontact shift field on the same grid as the spin 
%                density supplied.
%
% Note: minimal three-point schemes are used for the finite difference
%       operators. Increase stencil size if you have enough memory.
%
% Note: for further information on the equations and algorithms used in this
%       function see http://dx.doi.org/10.1039/C4CP03106G
%
% ilya.kuprov@weizmann.ac.uk
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=kpcs.m>

function [pcs_vals,pcs_cube]=kpcs(probden,chi,ranges,nxyz,method)

% Check consistency
grumble(probden,chi,ranges,nxyz);

% Isolate rank 2 component of chi
[~,~,rank2]=mat2sphten(chi);
chi=sphten2mat([],[],rank2);

% Compute axis extents
extents(1)=ranges(2)-ranges(1);
extents(2)=ranges(4)-ranges(3);
extents(3)=ranges(6)-ranges(5);

% Decide the method
switch method
    
    case 'fft'

        % Generate Fourier derivative multipliers
        [X,Y,Z]=ndgrid(fftdiff(1,size(probden,1),extents(1)/size(probden,1)),...
                       fftdiff(1,size(probden,2),extents(2)/size(probden,2)),...
                       fftdiff(1,size(probden,3),extents(3)/size(probden,3)));

        % Get Laplace operator in Fourier space
        L=X.^2+Y.^2+Z.^2; L(L==0)=1;

        % Get Kuprov operator in Fourier space
        K=-(1/3)*(chi(1,1)*X.*X+chi(1,2)*X.*Y+chi(1,3)*X.*Z+...
                  chi(2,1)*Y.*X+chi(2,2)*Y.*Y+chi(2,3)*Y.*Z+...
                  chi(3,1)*Z.*X+chi(3,2)*Z.*Y+chi(3,3)*Z.*Z);
              
        % Solve Kuprov equation (in ppm)
        pcs_cube=1e6*real(ifftn(K.*fftn(probden)./L));
       
    case 'fdiff'
        
        % Get 3-point Laplace operator
        L=fdlap(size(probden),extents,3);
        
        % Get 3-point Kuprov operator
        K=fdkup(size(probden),extents,chi,3);
        
        % Solve Kuprov equation (in ppm)
        pcs_cube=1e6*reshape(cgs(L,K*probden(:),1e-10,1000),size(probden));
        
    otherwise
        
        % Complain and bomb out
        error('unknown solver.');
        
end

% Deallocate variables
clear('L','K');

% Compute PCS values at the nuclei
[X,Y,Z]=ndgrid(linspace(ranges(1),ranges(2),size(pcs_cube,1)),...
               linspace(ranges(3),ranges(4),size(pcs_cube,2)),...
               linspace(ranges(5),ranges(6),size(pcs_cube,3)));
pcs_vals=interpn(X,Y,Z,pcs_cube,nxyz(:,1),nxyz(:,2),nxyz(:,3),'cubic');

end

% Consistency enforcement
function grumble(probden,chi,ranges,nxyz)
if (~isnumeric(nxyz))||(size(nxyz,2)~=3)||(~isreal(nxyz))
    error('nxyz parameter should be a real matrix with three columns.');
end
if (~isnumeric(chi))||(size(chi,1)~=3)||(size(chi,2)~=3)||(~isreal(chi))
    error('chi parameter should be a real 3x3 matrix.');
end
if (~isnumeric(probden))||(ndims(probden)~=3)||(~isreal(probden))
    error('probden parameter should be a real non-negative 3D array.');
end
if (~isnumeric(ranges))||(numel(ranges)~=6)||(~isreal(ranges))||...
   (ranges(1)>=ranges(2))||(ranges(3)>=ranges(4))||(ranges(5)>=ranges(6))
    error('ranges mult be a six-element array of reals with elements 1, 3, 5 smaller than 2, 4, 6 respectively.');
end
end

% The aim of science is to make difficult things understandable 
% in a simpler way; the aim of poetry is to state simple things 
% in an incomprehensible way. The two are incompatible.
%
% Paul Dirac

