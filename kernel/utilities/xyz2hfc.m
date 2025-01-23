% Converts point electron and nuclear coordinates into a hyper-
% fine interaction tensor. Syntax:
%
%                A=xyz2hfc(mxyz,exyz,isotope,nel)
%
% Parameters:
%
%     exyz     - Cartesian coordinates of the electron,
%                a three-element vector in Angstrom
%
%     nxyz     - Cartesian coordinates of the nucleus,
%                a three-element vector in Angstrom
%
%     isotope  - isitope specification, e.g. '13C'
%
%     nel      - number of unpaired electrons 
%
% Outputs:
%
%     A        - hyperfine coupling tensor, Gauss
%
% Note: Gauss units are used for hyperfine couplings because 
%       they do not depend on the electron g-tensor.
%
% ilya.kuprov@weizmann.ac.il
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=xyz2hfc.m>

function A=xyz2hfc(exyz,nxyz,isotope,nel)

% Check consistency
grumble(exyz,nxyz,isotope,nel);

% Fundamental constants
hbar=1.054571730e-34; 
mu0=4*pi*1e-7;

% Get magnetogyric ratios
gamma_n=spin(isotope);

% Set the origin
nxyz=nxyz-exyz;
        
% Collect fundamental constants
C=10^4*gamma_n*hbar*mu0*nel/(4*pi*(1e-10)^3);

% Compute the dipolar matrix
D=3*(nxyz'*nxyz)/norm(nxyz,2)^5-eye(3)/norm(nxyz,2)^3;
        
% Compute the dipolar coupling matrix
A=C*D;
        
end

% Consistency enforcement
function grumble(mxyz,nxyz,isotope,nel)
if (~isnumeric(mxyz))||(~isreal(mxyz))||(numel(mxyz)~=3)
    error('e_xyz must be a three-element real vector.');
end
if (~isnumeric(nxyz))||(~isreal(nxyz))||(numel(nxyz)~=3)
    error('n_xyz must be a three-element real vector.');
end
if ~all(size(mxyz)==size(nxyz))
    error('e_xyz and n_xyz must have the same dimension.');
end
if ~ischar(isotope)
    error('isotope specification must be a character string.');
end
if (~isnumeric(nel))||(~isreal(nel))||...
   (numel(nel)~=1)||(mod(nel,1)~=0)||(nel<1)
    error('nel must be a non-negative real integer.');
end
end

% "English people as a whole have a rooted distrust of total
%  abstainers as politicians."
%
% The Very Rev'd H. Hensley Henson, then Dean of Durham, 
% in a letter to the Times, 4th January 1916.

