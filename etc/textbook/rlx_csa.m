% Redfield theory expressions for CSA relaxation, including
% contributions from the antisymmetric part. Syntax:
%
%            [r1,r2]=rlx_csa(B0,isotope,Z,tau_c)
%
% Parameters:
%
%    B0       - magnet field, Tesla
%
%    isotope  - the spins involved, e.g. '15N'
%
%    Z        - chemical shift tensor, 3x3
%               matrix in ppm
%
%    tau_c    - second rank (1/6D) rotational 
%               correlation time, seconds
%
% Outputs:
%
%    r1       - longitudinal relaxation rate, Hz
%
%    r2       - transverse relaxation rate, Hz
%
% Note: CSA relaxation rate expressions do not depend on 
%       the spin quantum number
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=rlx_csa.m>

function [r1,r2]=rlx_csa(B0,isotope,Z,tau_c)

% Check consistency
grumble(B0,isotope,Z,tau_c);

% Get the carrier frequency
omega=B0*(1+trace(1e-6*Z)/3)*spin(isotope);

% Get Blicharski's invariants
[Lsq,Dsq]=blinv(1e-6*Z);

% Do the textbook math
r1=(1/2)*Lsq*omega^2*tau_c/(1+9*tau_c^2*omega^2)+...
   (2/15)*Dsq*omega^2*tau_c/(1+tau_c^2*omega^2);
r2=(1/4)*Lsq*omega^2*tau_c/(1+9*tau_c^2*omega^2)+...
   (1/45)*Dsq*omega^2*tau_c*(4+3/(1+tau_c^2*omega^2));

end

% Consistency enforcement
function grumble(B0,isotope,Z,tau_c)
if (~isnumeric(B0))||(~isreal(B0))||(~isscalar(B0))
    error('B0 must be a real number.');
end
if (~isnumeric(tau_c))||(~isreal(tau_c))||...
   (~isscalar(tau_c))||(tau_c<=0)
    error('tau_c must be a positive real number.');
end
if ~ischar(isotope)
    error('isotope must be a character string.');
end
if (~isnumeric(Z))||(~isreal(Z))||...
   (size(Z,1)~=3)||(size(Z,2)~=3)
    error('Z must be a real 3x3 matrix.');
end
end

% Reality is that which, when you stop believing 
% in it, doesn't go away.
%
% Philip K. Dick

