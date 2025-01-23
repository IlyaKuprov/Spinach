% Converts hyperfine coupling tensors and susceptibility tensors into
% pseudocontact shifts (contact component is not included) using Equa-
% tion 10 from http://dx.doi.org/10.1039/C4CP03106G. Syntax:
%
%              [pcs,pcs_tensor]=hfc2pcs(A,chi,isotope,nel)
%
% Parameters:
%
%            A    - hyperfine coupling tensor, Gauss
%
%          chi    - magnetic susceptibility tensor, Angstrom^3
%
%        isotope  - isotope i.e. '1H'
%
%           nel   - number of unpaired electrons 
%
% Outputs:
%
%          pcs    - isotropic pseudocontact shift, ppm
%
%   pcs_tensor    - pseudocontact shift tensor, ppm
%
% Note: Gauss units are used for hyperfine couplings because they do
%       not depend on the electron g-tensor.
%
% ilya.kuprov@weizmann.ac.il
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=hfc2pcs.m>

function [pcs,pcs_tensor]=hfc2pcs(A,chi,isotope,nel)

% Check consistency
grumble(A,chi,isotope,nel);

% Keep rank 2 to generate only the pseudocontact part
[~,~,rank2]=mat2sphten(A);     A=sphten2mat(0,[0 0 0],rank2);
[~,~,rank2]=mat2sphten(chi); chi=sphten2mat(0,[0 0 0],rank2);

% Fundamental constants
gamma_n=spin(isotope); 
hbar=1.05457173e-34; 
mu0=4*pi*1e-7;  

% Collect fundamental constants
C=10^4*gamma_n*hbar*mu0*nel/(4*pi*(1e-10)^3);

% Compute the full paramagnetic shift tensor
pcs_tensor=1e6*(1/(4*pi))*A*chi/C;

% Compute the isotropic part
pcs=trace(pcs_tensor)/3;

end

% Consistency enforcement
function grumble(A,chi,isotope,nel)
if (~isnumeric(A))||(~isreal(A))||...
   (~issymmetric(A))||(any(size(A)~=3))
    error('A must be a real symmetric 3x3 matrix.');
end
if (~isnumeric(chi))||(~isreal(chi))||...
   (~issymmetric(chi))||(any(size(chi)~=3))
    error('chi must be a real symmetric 3x3 matrix.');
end
if ~ischar(isotope)
    error('isotope must be a character string.');
end
if (~isnumeric(nel))||(~isreal(nel))||...
   (numel(nel)~=1)||(mod(nel,1)~=0)||(nel<1)
    error('nel must be a non-negative real integer.');
end
end

% The problem isn't that Johnny can't read. The problem isn't
% even that Johnny can't think. The problem is that Johnny 
% doesn't know what thinking is; he confuses it with feeling.
%
% Thomas Sowell

