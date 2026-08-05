% Coherent state of a bosonic mode. Builds the normalised trunca-
% tion of the coherent state with the specified amplitude on the
% specified bosonic mode, with unit operators on all other parti-
% cles of the system. Syntax:
%
%               rho=coherent(spin_system,mode,alpha)
%
% Parameters:
%
%    mode   - index of a bosonic mode in sys.isotopes
%
%    alpha  - coherent state amplitude, a complex scalar
%
% Outputs:
%
%    rho    - coherent state density matrix (zeeman-hilb)
%             or its vectorisation (zeeman-liouv)
%
% Note: the Fock space truncation of the mode chops the tail of
%       the Poisson distribution; the state is renormalised after
%       the truncation and the lost weight is reported.
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=coherent.m>

function rho=coherent(spin_system,mode,alpha)

% Check consistency
grumble(spin_system,mode,alpha);

% Fock state amplitudes up to the truncation level
nlev=spin_system.comp.mults(mode);
amps=alpha.^(0:(nlev-1))./sqrt(factorial(0:(nlev-1)));

% Report the Poisson distribution weight lost to the truncation
lost_weight=1-exp(-abs(alpha)^2)*sum(abs(amps).^2);
report(spin_system,['coherent state truncation lost ' ...
                    num2str(lost_weight,'%0.5g') ' of the norm.']);

% Normalise the truncated state
amps=amps/norm(amps,2);

% Build the mode density matrix
rho_mode=sparse(amps.'*conj(amps));

% Kron up the unit operators of the other particles
rho=1;
for n=1:spin_system.comp.nspins
    if n==mode
        rho=kron(rho,rho_mode);
    else
        rho=kron(rho,speye(spin_system.comp.mults(n)));
    end
end

% Convert into the current formalism
switch spin_system.bas.formalism

    case 'zeeman-hilb'

        % Density matrix is already built
        rho=full(rho);

    case 'zeeman-liouv'

        % Stretch the density matrix
        rho=rho(:);

    otherwise

        % Complain and bomb out
        error('coherent states are only available in Zeeman formalisms.');

end

end

% Consistency enforcement
function grumble(spin_system,mode,alpha)
if ~isfield(spin_system,'bas')
    error('basis set information is missing, run basis() before calling this function.');
end
if (~isnumeric(mode))||(~isscalar(mode))||(~isreal(mode))||...
   (mod(mode,1)~=0)||(mode<1)||(mode>spin_system.comp.nspins)
    error('mode must be the index of a particle in sys.isotopes.');
end
if ~ismember(spin_system.comp.types{mode},{'C','V','T'})
    error('the specified particle is not a bosonic mode.');
end
if (~isnumeric(alpha))||(~isscalar(alpha))||(~isfinite(alpha))
    error('alpha must be a finite complex scalar.');
end
end

