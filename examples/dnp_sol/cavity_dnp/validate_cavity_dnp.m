% Validate the full low-field E/209Bi/Cn DNP example in sphten-liouv.
% Compares Fock cutoffs, pulse Fokker-Planck ranks, RF saturation lengths,
% and reset timesteps. Independent native controls remove the RF pulse
% or cavity coupling, check lossless angular-momentum conservation, and
% compare bare cavity decay with the exact photon-loss exponential.
%
%                validation=validate_cavity_dnp(reference)
%
% reference is the result from cavity_dnp_minimal(3,2,300).
% A no-argument call runs that reference calculation first.
% No files or figures are written by this function.
%
% Calculation time: minutes; requests a 64-worker parallel pool
%
% ilya.kuprov@weizmann.ac.il

function validation=validate_cavity_dnp(reference)

% Reference calculation
if nargin==0
    reference=cavity_dnp_minimal(3,2,300);
end

% Spin system, reset generator, and acquisition parameters
spin_system=reference.spin_system;
H=reference.H;
R=reference.R;
observables=reference.observables;
rho=reference.rho;
time_step=reference.time(2);
nsteps=numel(reference.time)-1;

% Increase the photon cutoff without changing the physical model
cutoff=cavity_dnp_minimal(4,2,300);
cutoff_five=cavity_dnp_minimal(5,2,300);
cutoff_error=max([max(abs(cutoff.iz-reference.iz)) ...
                  max(abs(cutoff_five.iz-reference.iz)) ...
                  max(abs(cutoff_five.iz-cutoff.iz))]);
photon_cutoff_error=max([max(abs(cutoff.photons-reference.photons)) ...
                         max(abs(cutoff_five.photons-reference.photons))]);

% Increase the pulse phase-grid rank
phase_grid=cavity_dnp_minimal(3,3,300);
phase_grid_error=max(abs(phase_grid.iz-reference.iz));

% Double the physical saturation preparation length
saturation=cavity_dnp_minimal(3,2,600);
saturation_error=max([abs(saturation.initial_iz-reference.initial_iz) ...
                      max(abs(saturation.iz-reference.iz))]);

% Halve the reset time step
fine=evolution(spin_system,H+1i*R,observables,rho,...
               time_step/2,2*nsteps,'multichannel');
time_error=max(abs(real(fine(2,1:2:end))-reference.iz));

% Use a non-dyadic step ratio to test scaling-and-squaring errors
shorter_step=evolution(spin_system,H+1i*R,observables,rho,...
                       4*time_step/5,5*nsteps/4,'multichannel');
non_dyadic_error=max(abs(real(shorter_step(2,1:5:end))-reference.iz(1:4:end)));

% Omit the selective RF pulse
no_rf=evolution(spin_system,H+1i*R,observables,reference.rho_no_rf,...
                time_step,nsteps,'multichannel');
no_rf_gain=real(no_rf(2,end))-reference.initial_iz;

% Disconnect the cavity without changing photon loss
uncoupled_system=spin_system;
uncoupled_system.inter.modes.exchange{1,3}=0;
H_off=hamiltonian(uncoupled_system);
no_exchange=evolution(uncoupled_system,H_off+1i*R,observables,rho,...
                      time_step,nsteps,'multichannel');
no_exchange_photons=max(abs(no_exchange(3,:)));

% Compare one-photon decay with the exact exponential
unit=state(spin_system,'E',1);
one_photon=state(spin_system,'BL2',3);
one_photon=one_photon/(unit'*one_photon);
kappa=2*pi*210e3;
fast_time=linspace(0,10/kappa,101);
empty_decay=evolution(uncoupled_system,H_off+1i*R,observables,one_photon,...
                      fast_time(2),100,'multichannel');
photon_decay_error=max(abs(real(empty_decay(3,:))-exp(-kappa*fast_time)));

% Resolve the photon transient at two cavity cutoffs
photon_transient=evolution(spin_system,H+1i*R,observables,rho,...
                           fast_time(2),100,'multichannel');
fast_five=evolution(cutoff_five.spin_system,cutoff_five.H+1i*cutoff_five.R,...
                    cutoff_five.observables,cutoff_five.rho,...
                    fast_time(2),100,'multichannel');
fast_photon_error=max(abs(photon_transient(3,:)-fast_five(3,:)));
fast_photon_peak=max(real(photon_transient(3,:)));

% Native lossless conservation of total spin projection plus photon number
conserved=(operator(spin_system,'Lz',1,'left')+operator(spin_system,'Lz',2,'left')+...
           operator(spin_system,'N',3,'left'))'*unit;
lossless=evolution(spin_system,H,conserved,rho,1e-5,100,'observable');
conservation_error=max(abs(lossless-lossless(1)));

% Numerical convergence and physical controls
validation.cutoff_error=full(cutoff_error);
validation.photon_cutoff_error=full(photon_cutoff_error);
validation.phase_grid_error=full(phase_grid_error);
validation.saturation_error=full(saturation_error);
validation.time_error=full(time_error);
validation.non_dyadic_error=full(non_dyadic_error);
validation.fast_photon_error=full(fast_photon_error);
validation.fast_photon_peak=full(fast_photon_peak);
validation.no_rf_gain=full(no_rf_gain);
validation.no_exchange_photons=full(no_exchange_photons);
validation.photon_decay_error=full(photon_decay_error);
validation.conservation_error=full(conservation_error);
validation.trace_error=full(max(abs(reference.trace-1)));
validation.top_population=full(max(abs(reference.top_population)));
validation.initial_iz=reference.initial_iz;
validation.delta_iz=reference.delta_iz;
validation.scramble_change=reference.scramble_change;
validation.preparation_trace=reference.preparation_trace;
validation.transition_hz=reference.transition_hz;
validation.matlab_version=version;
disp(validation);

% Require convergence and the expected pumping mechanism
assert(strcmp(spin_system.bas.formalism,'sphten-liouv'));
assert(isequal(spin_system.comp.isotopes(1:2),{'E','209Bi'}));
assert(cutoff_error<1e-4&&photon_cutoff_error<1e-8);
assert(phase_grid_error<1e-4&&saturation_error<1e-4&&time_error<1e-4);
assert(non_dyadic_error<1e-4&&fast_photon_error<1e-9&&fast_photon_peak>1e-9);
assert(abs(no_rf_gain)<0.001&&reference.delta_iz>0.05);
assert(no_exchange_photons<1e-8&&photon_decay_error<1e-7);
assert(conservation_error<1e-6&&validation.trace_error<1e-5);
assert(validation.top_population<1e-8&&reference.scramble_change<1e-5);
assert(abs(reference.preparation_trace-1)<1e-5);
assert(min(reference.photons)>-1e-8&&max(abs(reference.iz))<=4.5);
fprintf('CAVITY_DNP_VALIDATION_SUCCESS\n');

end


