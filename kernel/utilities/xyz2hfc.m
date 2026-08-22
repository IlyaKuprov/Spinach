% Converts point electron and nuclear coordinates into a hyper-
% fine interaction tensor. Syntax:
%
%                  A=xyz2hfc(exyz,nxyz,isotope)
%
% Parameters:
%
%     exyz     - Cartesian coordinates of the electron,
%                a 1x3 row vector in Angstrom
%
%     nxyz     - Cartesian coordinates of the nucleus,
%                a 1x3 row vector in Angstrom
%
%     isotope  - isitope specification, e.g. '13C'
%
% Outputs:
%
%     A        - hyperfine coupling tensor, Gauss
%
% Note: Gauss units are used for hyperfine couplings because 
%       they do not depend on the electron g-tensor.
%
% Note: the tensor returned is the one that enters the spin
%       Hamiltonian as S*A*I; it does not scale with the num-
%       ber of unpaired electrons because the electron spin
%       operator already carries that magnitude.
%
% ilya.kuprov@weizmann.ac.il
% e.suturina@bath.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=xyz2hfc.m>

function A=xyz2hfc(exyz,nxyz,isotope)

% Check consistency
grumble(exyz,nxyz,isotope);

% Fundamental constants
hbar=1.054571730e-34; 
mu0=4*pi*1e-7;

% Get magnetogyric ratios
gamma_n=spin(isotope);

% Set the origin
nxyz=nxyz-exyz;
        
% Collect fundamental constants
C=10^4*gamma_n*hbar*mu0/(4*pi*(1e-10)^3);

% Compute the dipolar matrix
D=3*(nxyz'*nxyz)/norm(nxyz,2)^5-eye(3)/norm(nxyz,2)^3;
        
% Compute the dipolar coupling matrix
A=C*D;
        
end

% Consistency enforcement
function grumble(exyz,nxyz,isotope)
if (~isnumeric(exyz))||(~isreal(exyz))||(~isequal(size(exyz),[1 3]))
    error('exyz must be a 1x3 real row vector.');
end
if (~isnumeric(nxyz))||(~isreal(nxyz))||(~isequal(size(nxyz),[1 3]))
    error('nxyz must be a 1x3 real row vector.');
end
if norm(nxyz-exyz,2)==0
    error('e_xyz and n_xyz coordinates must be different.');
end
if ~ischar(isotope)
    error('isotope specification must be a character string.');
end
end

% "English people as a whole have a rooted distrust of total
%  abstainers as politicians."
%
% The Very Rev'd H. Hensley Henson, then Dean of Durham, 
% in a letter to the Times, 4th January 1916.

