% Analytical expression for a DEER trace for two spins in the 
% presence of dipolar and exchange coupling. Syntax:
%
%                    deer=deer_analyt(D,J,t)
%
% Parameters:
%
%     D - dipolar coupling, angular frequency units, the
%         coefficient in front of (1-3*cos(theta)^2)*Lz*Sz
%         in the spin Hamiltonian
%
%     J - exchange coupling, angular frequency units, NMR
%         convention (no factor of 2 in front), the coef- 
%         ficient in front of L*S in the spin Hamiltonian
%
%     t - array of time points, seconds
%
% Output:
%
%  deer - an array of DEER form factor values of the same
%         dimension as t
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=deer_analyt.m>

function deer=deer_analyt(D,J,t)

% Check consistency
grumble(D,J,t);

% Use Kuprov's formula (http://dx.doi.org/10.1038/ncomms14842)
deer=sqrt(pi./(6*D*t)).*(cos((D+J)*t).*fresnelc(sqrt(6*D*t/pi))+...
                         sin((D+J)*t).*fresnels(sqrt(6*D*t/pi)));
                     
% Remove zero time indeterminacy
deer(t==0)=1;

end

% Consistency enforcement
function grumble(D,J,t)
if (~isnumeric(D))||(~isreal(D))||(~isscalar(D))
    error('D must be a real number.');
end
if (~isnumeric(J))||(~isreal(J))||(~isscalar(J))
    error('J must be a real number.');
end
if (~isnumeric(t))||(~isreal(t))
    error('t must be an array of real numbers.');
end
end

% If nothing is certainly right, then of course it follows that
% nothing is certainly wrong.
%
% C.S. Lewis

