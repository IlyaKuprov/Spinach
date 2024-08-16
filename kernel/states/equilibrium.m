% Returns the thermal equilibrium state at the current temperature. If the
% anisotropic part and the orientation parameters are not given, uses the
% isotropic Hamiltonian, otherwise uses the full Hamiltonian at the speci-
% fied orientation. Syntax:
%
%               rho=equilibrium(spin_system,H,Q,euler_angles)
%
% Arguments:
%
%    H            -  isotropic part of the Hamiltonian left side pro-
%                    duct superoperator (in Liouville space) or Hamil-
%                    tonian (in Hilbert space).
%
%    Q            -  irreducible components of the anisotropic part
%                    of the Hamiltonian left side product superopera-
%                    tor (in Liouville space) or Hamiltonian (in Hil-
%                    bert space), as returned by hamiltonian.m
%
%    euler_angles -  a row vector of Euler angles (in radians) speci-
%                    fying the system orientation relative to the in-
%                    put orientation. If the angles are not supplied,
%                    only isotropic part of the Hamiltonian is used.
%
% Outputs:
%
%    rho          - thermal equilibrium density matrix (Hilbert spa-
%                   ce) or state vector (Liouville space).
%
% WARNING: Liouville space calculations must supply left side product su-
%          peroperators, not commutation superoperators.
%
% WARNING: assumptions supplied to the hamiltonian.m call that generates
%          H and Q must be 'labframe'.
%
% WARNING: spin system ground states are commonly degenerate; absolute
%          zero temperatures are not supported.
%
% ledwards@cbs.mpg.de
% i.kuprov@soton.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=equilibrium.m>

function rho=equilibrium(spin_system,H,Q,euler_angles)

% Account for the orientation
if nargin==4
    
    % Check consistency
    grumble(spin_system,H,Q,euler_angles);
    
    % Build anisotropic Hamiltonian
    H=H+orientation(Q,euler_angles);
    
elseif nargin==2
    
    % Check consistency
    grumble(spin_system,H);

elseif nargin==1

    % Build the isotropic Hamiltonian
    H=hamiltonian(assume(spin_system,'labframe'),'left');

    % Check consistency
    grumble(spin_system,H);
    
else
    
    % Complain and bomb out
    error('incorrect number of input arguments.');
    
end

% Get the temperature factor
beta_factor=spin_system.tols.hbar/(spin_system.tols.kbol*spin_system.rlx.temperature);

% Decide how to proceed
switch spin_system.bas.formalism
    
    % Liouville space formalisms
    case {'sphten-liouv','zeeman-liouv'}

        % Get the thermodynamic unit state
        switch spin_system.bas.formalism

            case 'sphten-liouv'

                % Unit population of T(0,0) state, normalisation is
                % such because prod(spin_system.comp.mults) can be-
                % come too large for double precision arithmetic
                unit=sparse(1,1,1,size(spin_system.bas.basis,1),1);

            case 'zeeman-liouv'

                % Stretched unit matrix, normalisation matched to
                % the Hilbert space because systems are small
                unit=speye(prod(spin_system.comp.mults)); unit=unit(:);

        end

        % Catch silly calls
        if norm(H*unit,1)<1e-10
            error('H and Q must be left side product superops, not commutation superops.');
        end

        % Propagate unit state in imaginary time
        rho=step(spin_system,H,unit,-1i*beta_factor);

        % Return to CPU if appropriate
        if isa(rho,'gpuArray'), rho=gather(rho); end

        % Catch failed calls
        if any(isnan(rho(:)))
            error('numerical accuracy issues, temperature too low - switch to Hilbert space.');
        end
        
        % Divide by partition function
        rho=rho/dot(unit,rho);
        
    % Hilbert space
    case {'zeeman-hilb'}
        
        % Estimate the norm
        mat_norm=abs(beta_factor)*cheap_norm(H);
            
        % Determine scaling and squaring parameters
        n_squarings=max([0 ceil(log2(mat_norm))]); scaling_factor=2^n_squarings;
            
        % Scale the temperature factor
        if scaling_factor>1, beta_factor=beta_factor/scaling_factor; end
            
        % Compute density matrix as a propagator in imaginary time
        rho=propagator(spin_system,H,-1i*beta_factor); rho=rho/trace(rho);
            
        % Move to GPU if appropriate
        if ismember('gpu',spin_system.sys.enable), rho=gpuArray(rho); end
            
        % Carefully square up
        if scaling_factor>1
            for n=1:n_squarings
                rho=rho*rho; rho=rho/trace(rho); % Numerical stability
                rho=clean_up(spin_system,rho,spin_system.tols.prop_chop);
            end
        end
            
        % Return to CPU if appropriate
        if isa(rho,'gpuArray'), rho=gather(rho); end
            
end
        
end
    
% Consistency enforcement
function grumble(spin_system,H,Q,euler_angles)
if ~isnumeric(H)
    error('Hamiltonian must be numeric.');
end
if isempty(spin_system.rlx.temperature)
    error('temperature must be specified in inter.temperature');
end
if spin_system.rlx.temperature==0
    error('absolute zero temperature is not supported.');
end
if nargin==4
    if ~iscell(Q)
        error('Q must be a cell array.');
    end
    if (~isnumeric(euler_angles))||(~isreal(euler_angles))||(numel(euler_angles)~=3)
        error('euler_angles must be a real three-element vector.');
    end
end
end

% Compassion is a wonderful thing. It's what one feels when one looks at
% a squashed caterpillar. An elevating experience. One can let oneself go
% and spread - you know, like taking a girdle off. You don't have to hold
% your stomach, your heart or your spirit up -- when you feel compassion. 
% All you have to do is look down. It's much easier. When you look up, you
% get a pain in the neck. [...] Oh, it has an antithesis - but such a hard,
% demanding one... admiration, Mrs. Jones, admiration. But that takes more
% than a girdle. So I say that anyone for whom we can't feel sorry is a
% vicious person.
%
% Ayn Rand, "The Fountainhead"

