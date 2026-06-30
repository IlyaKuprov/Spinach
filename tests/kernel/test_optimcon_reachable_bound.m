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
                       'reachable_bound() and oc_certify() must return rigorous one-spin speed-limit certificates.');

% Build a one-spin transfer problem
S=pauli(2);
drift=sparse(2,2);
controls={S.y};
rho_init=S.x;
rho_targ=S.z;
amplitude_bound=100;

% Too-short pulse: the remaining Bloch angle is pi/2-1
certificate=reachable_bound(drift,controls,rho_init,rho_targ,...
                            amplitude_bound,0.01);
result=test_close(result,'reachable time lower bound',certificate.time_lower_bound,...
                  pi/(2*amplitude_bound),1e-14,1e-14,...
                  'minimum time must follow from the Hamiltonian spectral width');
result=test_close(result,'reachable overlap upper bound',certificate.overlap_upper_bound,...
                  cos(pi/2-1),1e-14,1e-14,...
                  'overlap upper bound must reflect the unreachable residual angle');
result=test_true(result,'reachable fidelity bound below one',...
                 certificate.fidelity_upper_bound<1,...
                 'a pulse shorter than the speed-limit time must have a non-trivial bound');

% Long enough pulse: the certificate saturates at unity
certificate_long=reachable_bound(drift,controls,rho_init,rho_targ,...
                                 [-amplitude_bound amplitude_bound],0.02);
result=test_close(result,'reachable long overlap bound',certificate_long.overlap_upper_bound,...
                  1,1e-14,1e-14,...
                  'a long enough pulse must have a unit upper bound');

% Build a minimal control object and check the oc_certify wrapper
spin_system=local_spin_system();
spin_system.control.operators=controls;
spin_system.control.rho_init={rho_init};
spin_system.control.rho_targ={rho_targ};
spin_system.control.pwr_levels=amplitude_bound;
spin_system.control.pulse_dt=[0.005 0.005];
spin_system.control.drifts={{drift}};
certificate_oc=oc_certify(spin_system);
result=test_close(result,'oc_certify wrapper bound',certificate_oc.overlap_upper_bound,...
                  certificate.overlap_upper_bound,1e-14,1e-14,...
                  'oc_certify must pass an optimcon-style problem to reachable_bound');
result=test_true(result,'oc_certify empty comparison',...
                 isempty(certificate_oc.best_fidelity)&&isempty(certificate_oc.fidelity_gap),...
                 'oc_certify must leave GRAPE comparison fields empty without a waveform');

end


function spin_system=local_spin_system()

% Build a minimal quiet Spinach object for Hilbert-space helpers
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.sys.disable={};
spin_system.tols.liouv_zero=1e-14;
spin_system.tols.small_matrix=64;
spin_system.tols.dense_matrix=0.5;
spin_system.tols.prop_chop=1e-14;
spin_system.bas.formalism='zeeman-hilb';

end
