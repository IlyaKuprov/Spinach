% Stimulated echo diagnostics for the Mims ENDOR sequence. Syntax:
%
%    stim_echo=endor_mims_echo(spin_system,parameters,H,R,K)
%
% Parameters:
%
%     parameters.spins  - working spins, normally {'E'}; spe-
%                         cify multiplicity if electron spin
%                         is not 1/2, for example {'7E'} for
%                         gadolinium
%
%     parameters.electrons - a vector of integers specifying 
%                            which spins in sys.isotopes are
%                            electrons
%
%     parameters.tau       - the delay between the first two
%                            90-degree pulses of the Mims
%                            ENDOR sequence, seconds; 200e-9
%                            is typical
%
%     parameters.n_dur     - duration of the nuclear pulse,
%                            which is NOT ACTUALLY APPLIED
%                            within this pulse sequence,
%                            seconds; 50e-6 is typical
%
%     parameters.nsteps    - number of time steps to make
%                            in the detection period, which
%                            runs from 0 to 2*paramters.tau
%
%     H    - Hamiltonian matrix, received from the context
%            function, normally powder() in this case
%
%     R    - relaxation superoperator, received from the context
%            function, normally powder() in this case
%
%     K    - kinetics superoperator, received from the context
%            function, normally powder() in this case
%
% Outputs:
%
%     stim_echo - stimulated echo seen in the Mims ENDOR sequ-
%                 ence in the absence of the nuclear RF pulse
%
% ilya.kuprov@weizmann.ac.il
%
% <https://spindynamics.org/wiki/index.php?title=endor_mims_echo.m>

function stim_echo=endor_mims_echo(spin_system,parameters,H,R,K)

% Check consistency
grumble(spin_system,parameters,H,R,K);

% Compose Liouvillian
L=H+1i*R+1i*K;

% Ideal pulse operators on all electrons
Ex=operator(spin_system,'Lx',parameters.electrons);
Ey=operator(spin_system,'Ly',parameters.electrons);

% Ideal initial and detection states
rho0=state(spin_system,'Lz',parameters.electrons);
coil=state(spin_system,'L+',parameters.electrons);

% Ideal pi/2 pulse on the electrons
rho=step(spin_system,Ex,rho0,pi/2);

% Run stimulated echo delay
rho=step(spin_system,L,rho,parameters.tau);

% Ideal pi/2 pulse on the electrons
rho=step(spin_system,Ex,rho,pi/2);

% Delay corresponding to the missing nuclear pulse
rho=evolution(spin_system,L,[],rho,parameters.n_dur,1,'final');

% Ideal pi/2 pulse on the electrons
rho=step(spin_system,Ey,rho,-pi/2);

% Digitise the stimulated echo
stim_echo=evolution(spin_system,L,coil,rho,...
                    2*parameters.tau/parameters.nsteps,...
                    parameters.nsteps,'observable');

end

% Consistency enforcement
function grumble(spin_system,parameters,H,R,K)
if (~isnumeric(H))||(~isnumeric(R))||(~isnumeric(K))||...
   (~ismatrix(H))||(~ismatrix(R))||(~ismatrix(K))
    error('H, R and K arguments must be matrices.');
end
if (~all(size(H)==size(R)))||(~all(size(R)==size(K)))
    error('H, R and K matrices must have the same dimension.');
end
if ~ismember(spin_system.bas.formalism,{'sphten-liouv','zeeman-liouv'})
    error('this function is only available in Liouville space.');
end
if ~isfield(parameters,'n_dur')
    error('nuclear pulse duration must be specified in parameters.n_dur field.');
end
if (~isnumeric(parameters.n_dur))||(~isreal(parameters.n_dur))||...
   (~isscalar(parameters.n_dur))||(parameters.n_dur<=0)
    error('parameters.n_dur must be a positive real scalar.');
end
if ~isfield(parameters,'tau')
    error('echo delay must be specified in parameters.tau field.');
end
if (~isnumeric(parameters.tau))||(~isreal(parameters.tau))||...
   (~isscalar(parameters.tau))||(parameters.tau<=0)
    error('parameters.tau must be a positive real scalar.');
end
if ~isfield(parameters,'electrons')
    error('electrons must be enumerated in parameters.electrons field.');
end
if ~isfield(parameters,'nsteps')
    error('number of steps in the detection period must be specified in parameters.nsteps field.');
end
end

% "I see that you have made three spelling mistakes."
%
% Last words of Marquis de Favras after
% reading his death sentence before be-
% ing hanged in 1790.

