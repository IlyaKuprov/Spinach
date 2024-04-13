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
%    rho  - the state to be pumped, from state()
%
%    rate - pumping rate, Hz
%
% Outputs:
%
%    R    - modified relaxation superoperator
%
% Note: for the pumping to work correctly, the unit state population 
%       (first element) in the state vector that R will be acting on 
%        must be set to 1.
%
% Note: this function is only available in sphten-liouv formalism, and
%       may be called repeatedly if multiple states are pumped.
%
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=magpump.m>

function R=magpump(spin_system,R,rho,rate)

% Check consisttency
grumble(spin_system,R,rho,rate);

% Add pumping as a coupling to unit state
R(:,1)=R(:,1)+rate*rho;

end

% Consistency enforcement
function grumble(spin_system,R,rho,rate)
if (~isnumeric(R))||(~ismatrix(R))
    error('R must be a matrix.');
end
if (~isnumeric(rho))||(~iscolumn(rho))
    error('rho must be a column vector.');
end
if (~isnumeric(rate))||(~real(rate))||(~isscalar(rate))
    error('rate must be a real scalar.');
end
if ~ismember(spin_system.bas.formalism,{'sphten-liouv'})
    error('this function is only available in sphten-liouv formalism.');
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

