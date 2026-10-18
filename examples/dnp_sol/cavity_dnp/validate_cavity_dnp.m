% Validate native Spinach cavity DNP trajectories and independent controls.
% Run after adding Spinach and this folder to the MATLAB path. Returns
% numerical diagnostics; does not write files or open figures.
%
% ilya.kuprov@weizmann.ac.il

function validation=validate_cavity_dnp()

% Native cavity cutoffs and full-density consistency
results=cell(1,3);
for n=1:3
    results{n}=cavity_dnp_minimal(n+2);
    assert(results{n}.checks.trace_error<1e-8);
    assert(results{n}.checks.hermiticity_error<1e-8);
    assert(results{n}.checks.min_eigenvalue>-1e-8);
    assert(results{n}.checks.exact_error<1e-8);
    assert(results{n}.checks.max_top_population<1e-10);
    assert(strcmp(results{n}.spin_system.comp.isotopes{3},['C' num2str(n+2)]));
    assert(size(results{n}.final_density,1)==20*(n+2));
end
reference=results{1}; spin_system=reference.spin_system;

% Cavity cutoff convergence of both spin populations and photons
population_error=0; photon_error=0;
for n=2:3
    population_error=max(population_error,max(abs(...
        results{n}.populations-reference.populations),[],'all'));
    photon_error=max(photon_error,max(abs(...
        results{n}.photons-reference.photons)));
end
assert(population_error<1e-8); assert(photon_error<1e-8);

% Genuine half-step trajectory compared at every original grid point
step=reference.wait_time(2); nsteps=numel(reference.wait_time)-1;
rho=reference.initial_after_pulses;
fine=evolution(spin_system,reference.generator,[],rho,...
               step/2,2*nsteps,'trajectory');
step_error=max(abs(fine(:,1:2:end)-reference.wait_trajectory),[],'all');
assert(step_error<1e-8);

% No donor-cavity coupling means no pumping during the waiting stage
zero=evolution(spin_system,1i*reference.relaxation,[],rho,...
               reference.wait_time(end),1,'final');
p6=reference.projectors{1}; p15=reference.projectors{3};
zero_enhancement=9*real(p6(:)'*zero);
assert(abs(zero_enhancement-1)<1e-8);

% Lossless native evolution must conserve excitation and show Rabi exchange
g=reference.parameters.g; kappa=reference.parameters.kappa;
lossless_time=linspace(0,2*pi/g,81);
lossless=evolution(spin_system,g*reference.coupling_generator,[],rho,...
                   lossless_time(2),numel(lossless_time)-1,'trajectory');
number=reference.number; excitation=p15+number;
excitation_error=max(abs(real(excitation(:)'*lossless)-1/9));
rabi_error=max(abs(real(p15(:)'*lossless)-cos(g*lossless_time).^2/9));
assert(excitation_error<1e-8); assert(rabi_error<1e-8);

% Resolve the cavity transient without adiabatic cavity elimination
fast_time=linspace(0,12/kappa,101);
fast=evolution(spin_system,reference.generator,[],rho,...
               fast_time(2),numel(fast_time)-1,'trajectory');
zeta=sqrt((kappa/4)^2-g^2);
rslow=-g^2/(kappa/4+zeta); rfast=-kappa/4-zeta;
exact=(g/(2*zeta))^2*(exp(rslow*fast_time)-exp(rfast*fast_time)).^2/9;
fast_photon_error=max(abs(real(number(:)'*fast)-exact));
assert(fast_photon_error<1e-12);

% Report measured errors and native-system dimensions
validation=struct('population_cutoff_error',population_error,...
    'photon_cutoff_error',photon_error,'half_step_error',step_error,...
    'no_coupling_enhancement',zero_enhancement,...
    'lossless_excitation_error',excitation_error,'lossless_rabi_error',rabi_error,...
    'fast_photon_error',fast_photon_error,'enhancement',reference.enhancement,...
    'delta_iz',reference.delta_iz,'transition_hz',reference.transition_hz,...
    'checks',reference.checks,'fock_dimensions',[3 4 5],...
    'hilbert_dimensions',[60 80 100],'matlab_version',version);
disp(validation);
fprintf('CAVITY_DNP_VALIDATION_SUCCESS\n');

end


