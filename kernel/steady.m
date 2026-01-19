% Steady state under the repeated action by the same dissi-
% pative evolution propagator. Syntax:
%
%         rho=steady(spin_system,P,rho,tol,method)
%
% Parameters:
%
%    P - propagator, an exponential of the Liouvillian
%        that contains a thermalised relaxation super-
%        operator (inter.equilibrium='IME' or 'dibari')
%        or a product thereof (for example, from a re-
%        peating block of a pulse or a pulse sequence)
%
%    rho - optional initial guess for the steady state,
%          a good one can significantly accelerate this
%          function (leave empty otherwise); the first
%          element (it corresponds to the density mat-
%          rix trace in sphten-liouv) must always be 1
%
%    method - 'newton' (default) for the Newton-Raphson
%             steady state solver, 'squaring' for propa-
%             gator squaring (much more expensive, but
%             unconditionally numerically stable)
%
% Outputs:
%
%    rho - steady state under the repeated applicati-
%          on of the propagator P
%
% Note: only implemented for sphten-liouv formalism so far.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=steady.m>

function rho=steady(spin_system,P,rho,method)

% Default initial guess
if (~exist('rho','var'))||isempty(rho)
    rho=zeros([size(P,2) 1],'like',1i); rho(1)=1;
end

% Default method
if (~exist('method','var'))||isempty(method)
    method='newton';
end

% Check consistency
grumble(spin_system,P,rho);

% Pick the method
switch method

    case 'squaring'

        % Iteration stats
        cheap_diff=1; n_iter=1;

        % Keep going
        while cheap_diff>spin_system.tols.stst_tol

            % Compute the square
            Psq=clean_up(spin_system,P*P,spin_system.tols.prop_chop);

            % Compute the difference and close the loop
            cheap_diff=max(abs(Psq-P),[],'all'); 
            P=Psq; n_iter=n_iter+1;

            % Detect algorithm stagnation
            if n_iter>30, error('steady state convergence failure.'); end

        end
        
        % Compute the state
        rho=P*rho;

    case 'newton'

        % Get the Jacobian
        J=P-speye(size(P)); du=1;

        % Pre-factor the Jacobian
        [LF,UF,RP]=lu(J(2:end,2:end));

        % Iteration counter
        n_iter=0;

        % Newton iteration
        while norm(du,2)>spin_system.tols.stst_tol

            % Compute the residual
            r=P*rho-rho; r=r(2:end);

            % Re-use LU factors
            du=-UF\(LF\(RP*r));

            % Update the steady state
            rho(2:end)=rho(2:end)+du; n_iter=n_iter+1;

            % Detect algorithm stagnation
            if n_iter>10, error('steady state convergence failure.'); end

        end

    otherwise

        % Complain and bomb out
        error('unknown steady state calculation method.');

end

end

% Consistency enforcement
function grumble(spin_system,P,rho)
if ~strcmp(spin_system.bas.formalism,'sphten-liouv')
    error('steady state is only available for sphten-liouv formalism.');
end
if (~isnumeric(rho))||(~isnumeric(P))
    error('P and rho must be numeric.');
end
if size(P,1)~=size(P,2)
    error('P must be a square matrix.');
end
if (P(1,1)~=1)||(norm(P(1,2:end),2)~=0)
    error('P(1,:) must be [1 0 0 0 ...]');
end
if norm(P(2:end,1),2)==0
    error('the relaxation superoperator must be thermalised.');
end
if ~iscolumn(rho)
    error('rho must be a column vector.');
end
if rho(1)~=1
    error('rho(1) must be equal to 1.');
end
end

% People in Seville protested against mosquitoes
% spreading the West Nile fever.
%
% News section 
% of The Spectator

