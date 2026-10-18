% Validate the full low-field E/209Bi/Cn DNP example in sphten-liouv.
% Compares Fock cutoffs, pulse Fokker-Planck ranks, RF saturation lengths,
% and reset timesteps. Independent native controls remove the RF pulse
% or cavity coupling, check lossless angular-momentum conservation, and
% compare bare cavity decay with the exact photon-loss exponential.
%
%                validation=validate_cavity_dnp(reference)
%
% reference is optional; if omitted, the default example is executed.
% No files or figures are written by this function.
%
% Calculation time: minutes; requests a 64-worker parallel pool
%
% ilya.kuprov@weizmann.ac.il

function validation=validate_cavity_dnp(reference)

% Default full-donor example
if nargin==0, reference=cavity_dnp_minimal(3,2,300); end
s=reference.spin_system; H=reference.H; R=reference.R;
obs=reference.observables; rho=reference.rho;
dt=reference.time(2); nsteps=numel(reference.time)-1;

% Increase the photon cutoff without changing the physical model
cutoff=cavity_dnp_minimal(4,2,300);
cutoff_five=cavity_dnp_minimal(5,2,300);
cutoff_error=max([max(abs(cutoff.iz-reference.iz)) ...
                  max(abs(cutoff_five.iz-reference.iz)) ...
                  max(abs(cutoff_five.iz-cutoff.iz))]);
photon_cutoff_error=max([max(abs(cutoff.photons-reference.photons)) ...
                         max(abs(cutoff_five.photons-reference.photons))]);

% Increase the pulse phase-grid rank, retaining the same laboratory fields
phase_grid=cavity_dnp_minimal(3,3,300);
phase_grid_error=max(abs(phase_grid.iz-reference.iz));

% Double the physical saturation preparation length
saturation=cavity_dnp_minimal(3,2,600);
saturation_error=max([abs(saturation.initial_iz-reference.initial_iz) ...
                      max(abs(saturation.iz-reference.iz))]);

% Halve the laboratory-frame reset timestep with native multichannel output
fine=evolution(s,H+1i*R,obs,rho,dt/2,2*nsteps,'multichannel');
time_error=max(abs(real(fine(2,1:2:end))-reference.iz));

% A non-dyadic step ratio independently tests scaling-and-squaring effects
other_grid=evolution(s,H+1i*R,obs,rho,4*dt/5,5*nsteps/4,'multichannel');
non_dyadic_error=max(abs(real(other_grid(2,1:5:end))-reference.iz(1:4:end)));

% Removing only the RF pulse must suppress the persistent nuclear transfer
no_rf=evolution(s,H+1i*R,obs,reference.rho_no_rf,dt,nsteps,'multichannel');
no_rf_gain=real(no_rf(2,end))-reference.initial_iz;

% Remove exchange through the native interaction record and rebuild
s_off=s; s_off.inter.modes.exchange{1,3}=0;
H_off=hamiltonian(s_off);
no_exchange=evolution(s_off,H_off+1i*R,obs,rho,dt,nsteps,'multichannel');
no_exchange_photons=max(abs(no_exchange(3,:)));

% Check the native linewidth and Fock-state normalization independently
unit=state(s,'E',1); one_photon=state(s,'BL2',3);
one_photon=one_photon/(unit'*one_photon);
kappa=2*pi*210e3; fast_time=linspace(0,10/kappa,101);
empty_decay=evolution(s_off,H_off+1i*R,obs,one_photon,...
                      fast_time(2),100,'multichannel');
photon_decay_error=max(abs(real(empty_decay(3,:))-exp(-kappa*fast_time)));

% Resolve actual donor-generated photon build-up at two cavity cutoffs
fast=evolution(s,H+1i*R,obs,rho,fast_time(2),100,'multichannel');
fast_five=evolution(cutoff_five.spin_system,cutoff_five.H+1i*cutoff_five.R,...
    cutoff_five.observables,cutoff_five.rho,fast_time(2),100,'multichannel');
fast_photon_error=max(abs(fast(3,:)-fast_five(3,:)));
fast_photon_peak=max(real(fast(3,:)));

% Native lossless conservation of total spin projection plus photon number
conserved=(operator(s,'Lz',1,'left')+operator(s,'Lz',2,'left')+...
           operator(s,'N',3,'left'))'*unit;
lossless=evolution(s,H,conserved,rho,1e-5,100,'observable');
conservation_error=max(abs(lossless-lossless(1)));

% Summarize measured numerical and physical controls before enforcing gates
validation=struct('cutoff_error',full(cutoff_error),...
    'photon_cutoff_error',full(photon_cutoff_error),...
    'phase_grid_error',full(phase_grid_error),...
    'saturation_error',full(saturation_error),'time_error',full(time_error),...
    'non_dyadic_error',full(non_dyadic_error),...
    'fast_photon_error',full(fast_photon_error),'fast_photon_peak',full(fast_photon_peak),...
    'no_rf_gain',full(no_rf_gain),'no_exchange_photons',full(no_exchange_photons),...
    'photon_decay_error',full(photon_decay_error),...
    'conservation_error',full(conservation_error),...
    'trace_error',full(max(abs(reference.trace-1))),...
    'top_population',full(max(abs(reference.top_population))),...
    'initial_iz',reference.initial_iz,'delta_iz',reference.delta_iz,...
    'scramble_change',reference.scramble_change,...
    'preparation_trace',reference.preparation_trace,...
    'transition_hz',reference.transition_hz,'matlab_version',version);
disp(validation);

% Fail closed on formalism, convergence, and the expected pumping mechanism
assert(strcmp(s.bas.formalism,'sphten-liouv'));
assert(isequal(s.comp.isotopes(1:2),{'E','209Bi'}));
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

