% One-spin GRAPE certification example. A short bounded waveform is
% compared with a closed-form reachability bound for the same two-level
% control problem.
%
% Calculation time: seconds.
%
% ilya.kuprov@weizmann.ac.il

function oc_certify_one_spin_grape()

% Minimal quiet Spinach object for Hilbert-space optimal control
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.sys.disable={};
spin_system.tols.liouv_zero=1e-14;
spin_system.tols.small_matrix=64;
spin_system.tols.dense_matrix=0.5;
spin_system.tols.prop_chop=1e-14;
spin_system.bas.formalism='zeeman-hilb';

% Spin-half operators
S=pauli(2);

% Optimal-control setup
control.operators={S.y};
control.rho_init={S.x};
control.rho_targ={S.z};
control.pwr_levels=100;
control.pulse_dt=0.005*ones(1,2);
control.drifts={{sparse(2,2)}};
control.method='lbfgs';
control.max_iter=0;
control.penalties={'none'};
control.p_weights=0;
control.l_bound=-100;
control.u_bound=100;
control.plotting={};
spin_system=optimcon(spin_system,control);

% Candidate waveform and certificate
problem.waveform=[70 70];
certificate=oc_certify(spin_system,problem);

% Report the result
disp(certificate);

end

