% Computes pseudocontact shift from a point electron centre at the 
% nuclear coordinates supplied. Syntax:
%
%                    pred_pcs=ppcs(nxyz,mxyz,chi)
%
% Parameters: 
%
%     chi     - magnetic susceptibility tensor in cubic Angstroms
%               as a 3x3 matrix, or its five unique components 
%               ordered as
%                 
%                [chi(1,1) chi(1,2) chi(1,3) chi(2,2) chi(2,3)]
%
%     nxyz    - nuclear coordinates as [x y z] with multiple rows,
%               at which PCS is to be evaluated, in Angstroms.
%
%     sxyz    - susceptibility centre coordinates as [x y z], in 
%               Angstroms.
%
% Output:
% 
%     pcs     - predicted pseudocontact shift (in ppm) at each of 
%               the nuclei.
%
% ilya.kuprov@weizmann.ac.uk
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=ppcs.m>

function pcs=ppcs(nxyz,sxyz,chi)

% Form susceptibility tensor
if numel(chi)==5
    chi=[chi(1)     chi(2)     chi(3); 
         chi(2)     chi(4)     chi(5);
         chi(3)     chi(5) -chi(1)-chi(4)];
end

% Check consistency
grumble(nxyz,sxyz,chi);

% Put the metal at the origin
nxyz=nxyz-kron(sxyz,ones(size(nxyz,1),1));

% Convert coordinates to spherical
[r,theta,phi]=xyz2sph(nxyz(:,1),nxyz(:,2),nxyz(:,3));

% Get the irreducible components
[~,~,chi2]=qform2sph(chi);

% Preallocate the answer
pcs=zeros(numel(r),1);

% Sum the spherical harmonic expansion
for m=[2  1  0 -1 -2]
    pcs=pcs+(1/(4*pi))*chi2(-m+3)*(r.^(-3)).*spher_harmon(2,m,theta,phi);
end
 
% Convert to ppm and clean up
pcs=1e6*real(pcs);
 
end

% Consistency enforcement
function grumble(nxyz,sxyz,chi)
if (~isnumeric(nxyz))||(size(nxyz,2)~=3)||(~isreal(nxyz))
    error('nxyz must be a real matrix with three columns.');
end
if (~isnumeric(sxyz))||(size(sxyz,2)~=3)||(size(sxyz,1)~=1)||(~isreal(sxyz))
    error('sxyz must be a real row vector with three elements.');
end
if (~isnumeric(chi))||(size(chi,1)~=3)||(size(chi,2)~=3)||(~isreal(chi))
    error('chi must be a real 3x3 matrix or a vector of five real numbers');
end
end

% It's not true I had nothing on. I had the radio on.
%
% Marilyn Monroe

