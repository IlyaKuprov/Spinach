% Numerical controls for the exclusively sphten-liouv cavity DNP example.
% Checks native state normalization, pulse population swaps, Fock-cutoff
% and time-grid convergence, cavity-free and lossless controls, excitation
% conservation, and the resolved damped Jaynes-Cummings photon transient.
%
%                  validation=validate_cavity_dnp()
%
% Calculation time: seconds
%
% ilya.kuprov@weizmann.ac.il

function validation=validate_cavity_dnp()

% Cavity cutoff convergence of the delivered implementation
reference=cavity_dnp_minimal(3); truncation_error=0;
for nlevels=[4 5]
    candidate=cavity_dnp_minimal(nlevels);
    truncation_error=max(truncation_error,...
        max(abs(candidate.populations-reference.populations),[],'all'));
end

% Use native objects and observables without changing formalism
s=reference.spin_system; H=reference.H; R=reference.R;
rho=reference.rho; obs=reference.observables;
dt=reference.time(2); nsteps=numel(reference.time)-1;

% Half-step convergence with the native evolution routine
fine=evolution(s,H+1i*R,[],rho,dt/2,2*nsteps,'trajectory');
time_error=max(abs(real(obs'*fine(:,1:2:end))-...
    [reference.populations;reference.iz;reference.photons;...
     reference.top_population;reference.trace]),[],'all');

% Remove only the exchange interaction through its native specification
s_off=s; s_off.inter.modes.exchange{1,3}=0;
H_off=hamiltonian(s_off);
off=evolution(s_off,H_off+1i*R,[],rho,dt,nsteps,'trajectory');
off_error=max(abs(real(obs(:,1:3)'*off)-reference.after_pulses(1:3)),[],'all');

% Remove cavity damping to obtain coherent exchange
lossless=evolution(s,H,[],rho,dt,nsteps,'trajectory');
unit=state(s,'E','all');
excitation=(operator(s,'ZL2',1,'left')+operator(s,'N',3,'left'))'*unit;
excitation_trace=real(excitation'*lossless);
conservation_error=max(abs(excitation_trace-excitation_trace(1)));
g=2*pi*100; kappa=2*pi*2e3;
lossless_exact=0.5*cos(g*reference.time).^2;
lossless_error=max(abs(real(obs(:,3)'*lossless)-lossless_exact));

% Resolve the fast photon transient independently of the output grid
fast_time=linspace(0,10/kappa,201);
fast=evolution(s,H+1i*R,[],rho,fast_time(2),200,'trajectory');
zeta=sqrt((kappa/4)^2-g^2);
rslow=-g^2/(kappa/4+zeta); rfast=-kappa/4-zeta;
photon_exact=0.5*(g/(2*zeta))^2*...
             (exp(rslow*fast_time)-exp(rfast*fast_time)).^2;
photon_error=max(abs(real(obs(:,5)'*fast)-photon_exact));

% Physical and numerical acceptance gates
trace_error=max(abs(reference.trace-1));
pulse_error=max(abs(reference.after_pulses(1:3)-[0.5;0;0.5]));

% Compact machine-readable validation record
validation=struct('truncation_error',truncation_error,...
    'time_error',full(time_error),'off_error',full(off_error),...
    'conservation_error',full(conservation_error),'lossless_error',full(lossless_error),...
    'photon_error',full(photon_error),'trace_error',trace_error,...
    'pulse_error',pulse_error,'enhancement',reference.enhancement,...
    'delta_iz',reference.delta_iz,'fock_dimensions',[3 4 5],...
    'formalism',s.bas.formalism,'matlab_version',version);
disp(validation);

% Fail closed on all physical and numerical checks
assert(strcmp(s.bas.formalism,'sphten-liouv'));
assert(truncation_error<1e-8&&time_error<1e-8);
assert(off_error<1e-8&&conservation_error<1e-8);
assert(lossless_error<1e-7&&photon_error<1e-8);
assert(trace_error<1e-8&&pulse_error<1e-8);
assert(min(reference.populations,[],'all')>-1e-8);
assert(max(abs(reference.top_population))<1e-8);
assert(reference.enhancement>1.99&&reference.delta_iz>0.49);

fprintf('CAVITY_DNP_VALIDATION_SUCCESS\n');

end

