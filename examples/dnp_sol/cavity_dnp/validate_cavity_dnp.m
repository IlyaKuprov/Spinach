% Validate cavity truncation, density matrices, analytic decay, and controls.
% Run after adding Spinach and this folder to the MATLAB path.
% Returns validation metrics without writing files or opening figures.
%
% ilya.kuprov@weizmann.ac.il

function validation=validate_cavity_dnp()

% Run identical protocols with increasing photon cutoffs
results=cell(1,3);
for n=1:3
    results{n}=cavity_dnp_minimal(n+1);
    fprintf('CUTOFF %d enhancement %.12g trace %.3g exact %.3g\n',...
        n+1,results{n}.enhancement,results{n}.checks.trace_error,...
        results{n}.checks.exact_error);
end
reference=results{3};

% Convergence and physical invariants
population_error=0; photon_error=0;
for n=1:3
    population_error=max(population_error,max(abs(...
        results{n}.populations-reference.populations),[],'all'));
    photon_error=max(photon_error,max(abs(results{n}.photons-reference.photons)));
    assert(results{n}.checks.trace_error<1e-8);
    assert(results{n}.checks.hermiticity_error<1e-10);
    assert(results{n}.checks.min_eigenvalue>-1e-10);
    assert(results{n}.checks.exact_error<1e-8);
end
assert(population_error<1e-8);
assert(photon_error<1e-8);
assert(results{3}.checks.max_top_population<1e-12);
assert(abs(reference.enhancement-(2-exp(-8)))<1e-6);

% Time-step refinement against a single full-interval Spinach propagator
P=propagator(reference.spin_system,reference.generator,reference.wait_time(end));
final=P*reference.initial_after_pulses;
step_error=norm(final-reference.final_density(:));
assert(step_error<1e-8);

% No coupling control: after the ideal pulses, no DNP gain is possible
r=reference; dimension=size(r.final_density,1); cavity=weyl(dimension/3);
a=kron(speye(3),cavity.a); number=a'*a; identity=speye(dimension);
R=r.parameters.kappa*(kron(conj(a),a)-0.5*kron(identity,number)-...
    0.5*kron(number.',identity));
P=propagator(r.spin_system,1i*R,r.wait_time(end));
no_coupling=reshape(P*r.initial_after_pulses,dimension,dimension);
projector=kron(sparse(1,1,1,3,3),cavity.u);
no_coupling_enhancement=9*real(trace(projector*no_coupling));
assert(abs(no_coupling_enhancement-1)<1e-10);

% No loss control: one resonant vacuum Rabi oscillation conserves excitation
lower=sparse(1,3,1,3,3);
H=r.parameters.g*(kron(lower,cavity.c)+kron(lower',cavity.a));
L=kron(identity,H)-kron(H.',identity);
period=pi/r.parameters.g;
P=propagator(r.spin_system,L,period/100);
rho=r.initial_after_pulses;
excitation=kron(sparse(3,3,1,3,3),cavity.u)+number;
excitation_error=0; rabi_error=0;
for n=0:100
    density=reshape(rho,dimension,dimension);
    excitation_error=max(excitation_error,abs(trace(excitation*density)-1/9));
    excited=real(trace(kron(sparse(3,3,1,3,3),cavity.u)*density));
    rabi_error=max(rabi_error,abs(excited-cos(pi*n/100)^2/9));
    rho=P*rho;
end
assert(excitation_error<1e-10);
assert(rabi_error<1e-10);

% Resolve the fast cavity build-up independently of the slow plot grid
fast_time=linspace(0,20/reference.parameters.kappa,101);
P=propagator(reference.spin_system,reference.generator,fast_time(2));
rho=reference.initial_after_pulses;
kappa=reference.parameters.kappa; g=reference.parameters.g;
zeta=sqrt((kappa/4)^2-g^2);
rslow=-g^2/(kappa/4+zeta); rfast=-kappa/4-zeta;
fast_photon_error=0;
for n=1:numel(fast_time)
    density=reshape(rho,dimension,dimension);
    exact=(g/(2*zeta))^2*(exp(rslow*fast_time(n))-...
                               exp(rfast*fast_time(n)))^2/9;
    fast_photon_error=max(fast_photon_error,abs(trace(number*density)-exact));
    rho=P*rho;
end
assert(fast_photon_error<1e-12);

% Persist all numerical data and machine-readable validation
validation=struct('population_cutoff_error',population_error,...
    'photon_cutoff_error',photon_error,'step_error',step_error,...
    'no_coupling_enhancement',no_coupling_enhancement,...
    'lossless_excitation_error',excitation_error,'lossless_rabi_error',rabi_error,...
    'fast_photon_error',fast_photon_error,...
    'enhancement',reference.enhancement,'delta_iz',reference.delta_iz,...
    'transition_hz',reference.transition_hz,'checks',reference.checks,...
    'matlab_version',version);
disp(validation);
fprintf('CAVITY_DNP_VALIDATION_SUCCESS\n');

end


