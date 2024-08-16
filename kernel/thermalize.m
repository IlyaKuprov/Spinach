% Modifies the relaxation superoperator to drive the system to the user-
% specified target state (inhomogeneous master equation formalism) or to
% the equilibrium state of the lab frame Hamiltonian at the temperature
% provided by the user (DiBari-Levitt formalism). Syntax:
%
%           R=thermalize(spin_system,R,HLSPS,T,rho_eq,method)
%
% Parameters:
%
%     R       - symmetric negative definite relaxation super-
%               operator that drives the system towards the
%               zero state vector; this may be obtained from
%               relaxation.m if inter.equilibrium is 'zero'
%
%     HLSPS   - lab frame Hamiltonian left side product super-
%               operator, available from hamiltonian.m (also
%               call orientation.m if necessary); this is not
%               required for IME formalism (pass empty array)
%
%     T       - absolute temperature, not required for the 
%               IME formalism (pass empty array)
%
%     rho_eq  - thermal equilibrium state, not required for 
%               the DiBari-Levitt formalism (pass empty array)
%
%     method  - 'dibari' for DiBari-Levitt thermalisation,
%               'IME' for the inhomogeneous master equation
%
% Outputs:
%
%     R       - thermalized relaxation superoperator
%
% Note: to work correctly, IME requires the population of the unit state 
%       in the state vector to be exactly 1. Spinach has no way of check-
%       ing or enforcing this requirement - take due care.
%
% Note: DiBari-Levitt method is computationally expensive, but tends to
%       work better than IME, particularly in exotic regimes.
%
% i.kuprov@soton.ac.uk
% fije@inano.au.dk
%
% <https://spindynamics.org/wiki/index.php?title=thermalize.m>

function R=thermalize(spin_system,R,HLSPS,T,rho_eq,method)

% Check consistency
grumble(spin_system,R,HLSPS,T,rho_eq,method);

% Choose the method
switch method
    
    case 'IME'

        % This is formalism-dependent
        switch spin_system.bas.formalism

            case 'sphten-liouv'

                % Unit state has unit population of T(0,0) state
                U=sparse(1,1,1,size(spin_system.bas.basis,1),1);

            case 'zeeman-liouv'

                % Unit state is a stretched unit matrix
                U=speye(prod(spin_system.comp.mults)); U=U(:);

            otherwise

                % Complain and bomb out
                error('this function is only available in Liouville space.');

        end
        
        % Apply IME correction
        R=R-kron(U',R*rho_eq);
        
    case 'dibari'
        
        % Get the temperature factor
        beta=spin_system.tols.hbar/(spin_system.tols.kbol*T);

        % Modify the relaxation superoperator
        R=R*propagator(spin_system,HLSPS,1i*beta);
        
    otherwise
        
        % Complain and bomb out
        error('unknown thermalization method.');
        
end

end

% Consistency enforcement
function grumble(spin_system,R,HLSPS,T,rho_eq,method)
if (~isnumeric(R))||(size(R,1)~=size(R,2))
    error('R must be a square matrix.');
end
unit=unit_state(spin_system);
if norm(R*unit,2)>1e-10
    error('R appears to be thermalized already.');
end
if ~ischar(method)
    error('method must be a character string.');
end
if strcmp(method,'IME')
    if isempty(rho_eq)
        error('rho_eq cannot be empty for IME formalism.');
    end
    if (~isnumeric(rho_eq))||(~iscolumn(rho_eq))
        error('rho_eq must be a column vector.');
    end
end
if strcmp(method,'dibari')
    if isempty(HLSPS)
        error('HLSPS cannot be empty for DiBari-Levitt formalism.');
    end
    if (~isnumeric(HLSPS))||(size(HLSPS,1)~=size(HLSPS,2))
        error('HLSPS must be a square matrix.');
    end
    if norm(HLSPS*unit,2)<1e-8
        error('HLSPS appears to be a commutation superoperator.');
    end
    if isempty(T)
        error('T cannot be empty for DiBari-Levitt formalism.');
    end
    if (~isnumeric(T))||(~isreal(T))||(~isscalar(T))||(T<=0)
        error('T must be a positive real scalar.');
    end
end
end

% Tretyakov, the owner of the famous picture gallery in St Petersburg, had
% ordered the guards not to let Ilya Repin, a famous painter, into the gal-
% lery after Repin was repeatedly seen turning up with brushes and paints,
% and making small fixes to his works that the gallery had bought.

