% Tests optimal-control reachability certificates. Syntax:
%
%                    result=test_optimcon_reachable_bound()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% ilya.kuprov@weizmann.ac.il

function result=test_optimcon_reachable_bound()

% Announce the test target
fprintf('TESTING: Optimal-control reachability certificate\n');

% State the test target
result=new_test_result('optimcon/reachable_bound',...
                       'Optimal-control reachability certificate',...
                       'reachable_bound() and oc_certify() must return rigorous Hilbert-space speed-limit certificates.');

% Build a two-spin Hilbert-space transfer problem
S=pauli(2);
E=speye(2);
drift=sparse(4,4);
controls={kron(S.y,E)};
rho_init=kron(S.x,E);
rho_targ=kron(S.z,E);
amplitude_bound=100;

% Too-short pulse: the remaining Hilbert-Schmidt angle is pi/2-1
certificate=reachable_bound(drift,controls,rho_init,rho_targ,...
                            amplitude_bound,0.01);
result=test_close(result,'reachable time lower bound',certificate.time_lower_bound,...
                  pi/(2*amplitude_bound),1e-14,1e-14,...
                  'minimum time must follow from the Hamiltonian spectral width');
result=test_close(result,'reachable overlap upper bound',certificate.overlap_upper_bound,...
                  cos(pi/2-1),1e-14,1e-14,...
                  'overlap upper bound must reflect the unreachable residual angle');
result=test_true(result,'reachable arbitrary dimension',...
                 strcmp(certificate.scope,'arbitrary-dimension Hilbert-space speed-limit certificate'),...
                 'the certificate must not be restricted to two-level systems');
result=test_true(result,'reachable fidelity bound below one',...
                 certificate.fidelity_upper_bound<1,...
                 'a pulse shorter than the speed-limit time must have a non-trivial bound');
result=test_true(result,'reachable SDP slots',...
                 isfield(certificate,'specialised_sdp_relaxations')&&...
                 isfield(certificate.specialised_sdp_relaxations,'single_excitation_chain'),...
                 'the certificate must reserve slots for specialised SDP relaxations');

% Long enough pulse: the certificate saturates at unity
certificate_long=reachable_bound(drift,controls,rho_init,rho_targ,...
                                 [-amplitude_bound amplitude_bound],0.02);
result=test_close(result,'reachable long overlap bound',certificate_long.overlap_upper_bound,...
                  1,1e-14,1e-14,...
                  'a long enough pulse must have a unit upper bound');

% Check Fubini-Study state-vector mode
psi_init=[1;0;0;0];
psi_targ=[0;0;1;0];
certificate_vec=reachable_bound(drift,controls,psi_init,psi_targ,...
                                amplitude_bound,0.01);
result=test_close(result,'state-vector speed bound',certificate_vec.time_lower_bound,...
                  pi/amplitude_bound,1e-14,1e-14,...
                  'state-vector certificates must use the Fubini-Study speed limit');

% Build a minimal control object and check the oc_certify wrapper
spin_system=local_spin_system([0.005 0.005]);
spin_system.control.operators=controls;
spin_system.control.rho_init={rho_init};
spin_system.control.rho_targ={rho_targ};
spin_system.control.pwr_levels=amplitude_bound;
spin_system.control.l_bound=-1;
spin_system.control.u_bound=1;
spin_system.control.drifts={{drift}};
certificate_oc=oc_certify(spin_system,struct());
result=test_close(result,'oc_certify power-scaled bound',certificate_oc.time_lower_bound,...
                  certificate.time_lower_bound,1e-14,1e-14,...
                  'oc_certify must scale optimcon bounds by pwr_levels');
result=test_true(result,'oc_certify empty comparison',...
                 isempty(certificate_oc.best_fidelity)&&isempty(certificate_oc.fidelity_gap),...
                 'oc_certify must leave waveform comparison fields empty without a waveform');

% A zero waveform must produce a finite certificate, not a GRAPE error
problem.waveform=zeros(1,2);
certificate_zero=oc_certify(spin_system,problem);
result=test_close(result,'zero waveform fidelity',certificate_zero.best_fidelity,...
                  0,1e-14,1e-14,...
                  'zero-overlap candidate waveforms must be accepted for comparison');
result=test_close(result,'zero waveform fidelity gap',certificate_zero.fidelity_gap,...
                  certificate_zero.fidelity_upper_bound,1e-14,1e-14,...
                  'fidelity gaps must use the same normalised scale');

% Time-dependent drift slices must all enter the conservative bound
spin_system.control.drifts={{drift 10*controls{1}}};
certificate_drift=oc_certify(spin_system,struct());
result=test_true(result,'time-dependent drift bound',...
                 certificate_drift.hamiltonian_spectral_width_bound>...
                 certificate_oc.hamiltonian_spectral_width_bound,...
                 'all time-dependent drift slices must contribute to the speed bound');

% Multi-objective controls must not be silently reduced to the first pair
spin_system.control.rho_init={rho_init,rho_init};
spin_system.control.rho_targ={rho_targ,rho_targ};
caught=false;
try
    oc_certify(spin_system,struct());
catch
    caught=true;
end
result=test_true(result,'multi-objective rejection',caught,...
                 'multi-objective controls must be rejected explicitly');

% Non-finite certificate inputs must be rejected
caught=false;
try
    reachable_bound(drift,controls,rho_init,rho_targ,Inf,0.01);
catch
    caught=true;
end
result=test_true(result,'non-finite input rejection',caught,...
                 'non-finite amplitude bounds must not produce certificates');

end


function spin_system=local_spin_system(pulse_dt)

% Build a minimal quiet Spinach object for Hilbert-space helpers
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.sys.disable={};
spin_system.tols.liouv_zero=1e-14;
spin_system.tols.small_matrix=64;
spin_system.tols.dense_matrix=0.5;
spin_system.tols.prop_chop=1e-14;
spin_system.bas.formalism='zeeman-hilb';

% Add minimal optimal-control metadata
spin_system.control.pulse_dt=pulse_dt;
spin_system.control.integrator='rectangle';
spin_system.control.dead_time=0;
spin_system.control.prefix=[];
spin_system.control.suffix=[];
spin_system.control.keyholes=cell(1,numel(pulse_dt));

end
