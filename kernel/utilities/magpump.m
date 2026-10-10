% Adds phenomenological pumping terms to the relaxation superoperator 
% to enable approximate simulation of CIDNP, PHIP and DNP type effects. 
% Syntax:
%
%                  R=magpump(spin_system,R,rho,rate)
%
% Parameters:
%
%    R    - relaxation superoperator, from relaxation()
%
%    rho  - unweighted polarisation shape to pump, from coil_state()
%
%    rate - pumping rate, Hz
%
% Outputs:
%
%    R    - modified relaxation superoperator
%
% Note: each substance is pumped through its own unit coordinate at
%       bas.offsets(n)+1. That population in the state vector on which
%       R acts supplies the instantaneous concentration, so pumping
%       scales with population without division by concentrations.
%
% Note: this function is only available in sphten-liouv formalism, and
%       may be called repeatedly if multiple states are pumped.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=magpump.m>

function R=magpump(spin_system,R,rho,rate)

% Check consistency
grumble(spin_system,R,rho,rate);

% Couple each local target to the unit coordinate of its own substance
for n=1:spin_system.bas.nsubst
    rows=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
    unit=spin_system.bas.offsets(n)+1;
    R(rows,unit)=R(rows,unit)+rate*rho(rows);
end

end

% Consistency enforcement
function grumble(spin_system,R,rho,rate)
if (~isnumeric(R))||(~ismatrix(R))
    error('R must be a matrix.');
end
if (~isnumeric(rho))||(~iscolumn(rho))
    error('rho must be a column vector.');
end
if (~isnumeric(rate))||(~isreal(rate))||(~isscalar(rate))||(~isfinite(rate))
    error('rate must be a finite real scalar.');
end
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available in sphten-liouv formalism.');
end
if any(rho(spin_system.bas.offsets(1:end-1)+1)~=0)
    error('unit state cannot be pumped.');
end
end

% I know of scarcely anything so apt to impress the imagination as the
% wonderful form of cosmic order expressed by the [Central Limit Theorem].
% The law would have been personified by the Greeks and deified, if they
% had known of it. It reigns with serenity and in complete self-effacement,
% amidst the wildest confusion. The huger the mob, and the greater the
% apparent anarchy, the more perfect is its sway. It is the supreme law of
% Unreason. Whenever a large sample of chaotic elements are taken in hand
% and marshalled in the order of their magnitude, an unsuspected and most
% beautiful form of regularity proves to have been latent all along.
% 
% Sir Francis Galton

