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
%     rho_eq  - unit-concentration target in each substance block,
%               from equilibrium with chem.concs set to ones; not required for
%               the DiBari-Levitt formalism (pass empty array)
%
%     method  - 'dibari' for DiBari-Levitt thermalisation,
%               'IME' for the inhomogeneous master equation
%
% Outputs:
%
%     R       - thermalized relaxation superoperator; in zeeman-hilb
%               IME, @(t,rho) returns the matrix derivative instead
%
% Note: IME is applied independently to each substance block. The target
%       states in this construction are unweighted, with unit population
%       at every substance unit coordinate. Acting on weighted states
%       scales the target by the instantaneous population, including zero.
%       Cross-substance blocks of R are not supported in IME.
%       Hilbert IME takes R in vectorised per-substance Liouville blocks
%       (dimension sum(D_n^2)) and a block-diagonal unit-trace rho_eq. Its
%       output is a matrix RHS, not a Hamiltonian or a step generator.
%
% Note: DiBari-Levitt method is computationally expensive, but tends to
%       work better than IME, particularly in exotic regimes.
%
% ilya.kuprov@weizmann.ac.il
% fije@inano.au.dk
%
% <https://spindynamics.org/wiki/index.php?title=thermalize.m>

function R=thermalize(spin_system,R,HLSPS,T,rho_eq,method)

% Check consistency
grumble(spin_system,R,HLSPS,T,rho_eq,method);

% Expose Hilbert IME as a matrix derivative rather than a Hamiltonian
if strcmp(spin_system.bas.formalism,'zeeman-hilb')&&strcmp(method,'IME')
    blocks=cell(spin_system.bas.nsubst,1);
    for n=1:spin_system.bas.nsubst
        idx=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
        blocks{n}=rho_eq(idx,idx);
    end
    zeeman=spin_system; zeeman.bas.formalism='zeeman-liouv';
    zeeman.bas.nstates=spin_system.bas.nstates.^2;
    zeeman.bas.offsets=[0;cumsum(zeeman.bas.nstates)];
    R=thermalize(zeeman,R,[],[],hilb2liouv(blocks,'statevec'),'IME');
    R=@(t,rho)hilb_action(spin_system,R,t,rho);
    return
end

% Choose the method
switch method
    
    case 'IME'

        % Apply IME independently within each substance block
        for n=1:spin_system.bas.nsubst
            idx=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
            switch spin_system.bas.formalism
                case 'sphten-liouv'
                    U=sparse(1,1,1,numel(idx),1);
                case 'zeeman-liouv'
                    U=speye(prod(spin_system.comp.mults(spin_system.chem.parts{n})));
                    U=U(:);
                otherwise
                    error('this function is only available in Liouville space.');
            end
            R(idx,idx)=R(idx,idx)-(R(idx,idx)*rho_eq(idx))*U';
        end

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
if strcmp(spin_system.bas.formalism,'zeeman-wavef')
    error('Spinach:thermalize:wavefunction',...
          'thermalisation is not supported in zeeman-wavef formalism.');
end
if (~isnumeric(R))||(size(R,1)~=size(R,2))
    error('R must be a square matrix.');
end
if strcmp(spin_system.bas.formalism,'zeeman-hilb')&&ischar(method)&&strcmp(method,'IME')
    dim=spin_system.bas.offsets(end);
    if size(R,1)~=sum(spin_system.bas.nstates.^2)||any(~isfinite(R),'all')
        error('Spinach:thermalize:hilbertMap','Hilbert IME requires a finite direct-sum Liouville relaxation map.');
    end
    if ~isnumeric(rho_eq)||~isequal(size(rho_eq),[dim dim])||any(~isfinite(rho_eq),'all')
        error('Spinach:thermalize:hilbertTarget','Hilbert IME requires a finite block-diagonal target matrix.');
    end
    for n=1:spin_system.bas.nsubst
        idx=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
        if nnz(rho_eq(idx,:))~=nnz(rho_eq(idx,idx))
            error('Spinach:thermalize:hilbertTarget','Hilbert IME requires a finite block-diagonal target matrix.');
        end
    end
    return
end
unit_system=spin_system; unit_system.chem.concs(:)=1;
unit=unit_state(unit_system);
if norm(R*unit,2)>1e-10
    error('R appears to be thermalized already.');
end
if ~ischar(method)
    error('method must be a character string.');
end
if ~ismember(method,{'IME','dibari'})
    error('method must be ''IME'' or ''dibari''.');
end
if strcmp(method,'IME')
    if spin_system.bas.nsubst>1
        for n=1:spin_system.bas.nsubst
            idx=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
            if nnz(R(idx,:))~=nnz(R(idx,idx))
                error('Spinach:thermalize:crossSubstanceRelaxation',...
                      'IME relaxation must not contain cross-substance blocks (substance %d).',n);
            end
        end
    end
    if isempty(rho_eq)
        error('rho_eq cannot be empty for IME formalism.');
    end
    if (~isnumeric(rho_eq))||(~iscolumn(rho_eq))
        error('rho_eq must be a column vector.');
    end
    if size(rho_eq,1)~=size(R,1)
        error('rho_eq and R must have matching dimensions.');
    end
    for n=1:spin_system.bas.nsubst
        idx=(spin_system.bas.offsets(n)+1):spin_system.bas.offsets(n+1);
        population=unit(idx)'*rho_eq(idx);
        if strcmp(spin_system.bas.formalism,'zeeman-liouv')
            population=sqrt(prod(spin_system.comp.mults(spin_system.chem.parts{n})))*population;
        end
        if abs(population-1)>1e-10
            error('Spinach:thermalize:targetConcentration',...
                  'IME requires unit-concentration target blocks; request equilibrium with chem.concs set to ones.');
        end
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

