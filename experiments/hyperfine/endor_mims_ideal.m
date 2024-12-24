% Mims ENDOR sequence with ideal electron pulses. Syntax:
%
%     endor_spec=endor_mims_ideal(spin_system,parameters,H,R,K)
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
%     endor_spec - Mims ENDOR spectrum, a vector of the same
%                  size as parameters.n_frq
%
% ilya.kuprov@weizmann.ac.uk
%
% <https://spindynamics.org/wiki/index.php?title=endor_mims_ideal.m>

function endor_spec=endor_mims_ideal(spin_system,parameters,H,R,K)

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

% Nutation frequencies at the specified RF B1 field
nut_freqs=-spin_system.inter.gammas*parameters.rf_b1_field;

% Nuclear pulse operators
Nx=sparse(0); Ny=sparse(0);
for n=parameters.nuclei

    % Gamma-weighted because isotopes may differ
    Nx=Nx+nut_freqs(n)*operator(spin_system,'Lx',n);
    Ny=Ny+nut_freqs(n)*operator(spin_system,'Ly',n);

end

% Ideal pi/2 pulse on the electrons
rho=step(spin_system,Ex,rho0,pi/2);

% Run stimulated echo delay
rho=step(spin_system,L,rho,parameters.tau);

% Ideal pi/2 pulse on the electrons
rho=step(spin_system,Ex,rho,pi/2);

% Kron up to make array over radiofrequencies
rho=kron(rho,ones(1,numel(parameters.n_frq)));

% Loop over the radiofrequency array
parfor n=1:numel(parameters.n_frq)                                  %#ok<*PFBNS>

    % Blast the nuclei and subtract the background
    rho(:,n)=shaped_pulse_af(spin_system,L,Nx,Ny,rho(:,n),parameters.n_frq(n),1, ...
                             parameters.n_dur,0,parameters.n_rnk,'expm')-        ...
             shaped_pulse_af(spin_system,L,Nx,Ny,rho(:,n),parameters.n_frq(n),0, ...
                             parameters.n_dur,0,parameters.n_rnk,'expm');       

end

% Ideal pi/2 pulse on the electrons
rho=step(spin_system,Ey,rho,-pi/2);

% Run stimulated echo delay
rho=step(spin_system,L,rho,parameters.tau);

% Detect echo tips
endor_spec=coil'*rho;
                                           
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
if ~isfield(parameters,'n_frq')
    error('nuclear pulse frequencies must be specified in parameters.n_frq field.');
end
if (~isnumeric(parameters.n_frq))||(~isreal(parameters.n_frq))
    error('parameters.n_frq must be an array of real numbers.');
end
if ~isfield(parameters,'rf_b1_field')
    error('RF B1 field must be specified in parameters.rf_b1_field variable.');
end
if (~isnumeric(parameters.rf_b1_field))||(~isreal(parameters.rf_b1_field))||...
   (~isscalar(parameters.rf_b1_field))
    error('parameters.rf_b1_field must be a real scalar.');
end
if ~isfield(parameters,'n_rnk')
    error('nuclear pulse grid rank must be specified in parameters.n_rnk field.');
end
if (~isnumeric(parameters.n_rnk))||(~isreal(parameters.n_rnk))||...
   (~isscalar(parameters.n_rnk))||(mod(parameters.n_rnk,1)~=0)||...
   (parameters.n_rnk<1)
    error('parameters.n_rnk must be a positive real integer.');
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
if ~isfield(parameters,'nuclei')
    error('nuclei must be enumerated in parameters.nuclei field.');
end
end

% "There is considerable overlap between the intelligence
%  of the smartest bears and the dumbest tourists."
%
% A forest ranger at the Yosemite national Park,
% about why it is hard to design the perfect gar-
% bage bin to keep bears from breaking into it.

