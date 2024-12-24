% Converts hyperfine coupling tensors and susceptibility tensors into
% paramagnetic shifts (contact + pseudocontact component) using Equa-
% tion 10 from http://dx.doi.org/10.1039/C4CP03106G. Syntax:
%
%             [pms,pms_tensor]=hfc2pms(A,chi,isotope,nel)
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
%          pms    - isotropic paramagnetic shift, ppm
%
%   pms_tensor    - paramagnetic shift tensor, ppm
%
% Note: Gauss units are used for hyperfine couplings because they do
%       not depend on the electron g-tensor.
%
% ilya.kuprov@weizmann.ac.uk
% e.suturina@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=hfc2pms.m>

function [pms,pms_tensor]=hfc2pms(A,chi,isotope,nel)

% Check consistency
grumble(A,chi,isotope,nel);

% Fundamental constants
gamma_n=spin(isotope); 
hbar=1.05457173e-34; 
mu0=4*pi*1e-7;  

% Collect fundamental constants
C=10^4*gamma_n*hbar*mu0*nel/(4*pi*(1e-10)^3);

% Compute the full paramagnetic shift tensor
pms_tensor=1e6*(1/(4*pi))*A*chi/C;

% Compute the isotropic part
pms=trace(pms_tensor)/3;

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

% If you tell a beautiful woman that she is beautiful, what have you given her? 
% It's no more than a fact and it has cost you nothing. But if you tell an ugly 
% woman that she is beautiful, you offer her the great homage of corrupting the 
% concept of beauty [...] you sacrifice your conscience, your reason, your inte-
% grity and your invaluable self-esteem.
%
% Ayn Rand, "Atlas Shrugged"

