% Tests one-spin optimal-control setup and Hilbert-space GRAPE. Syntax:
%
%                    result=test_optimcon_grape_one_spin()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that optimcon() accepts a minimal one-spin Hilbert-space
% control problem, and that grape_hilb() fidelity and gradient agree with
% independent matrix exponentiation and finite differences.
%
% ilya.kuprov@weizmann.ac.il

function result=test_optimcon_grape_one_spin()

% Announce the test target
fprintf('TESTING: One-spin optimal-control setup and GRAPE gradient\n');

% State the optimal-control target of the test
result=new_test_result('optimcon/grape_one_spin',...
                       'One-spin optimal-control setup and GRAPE gradient',...
                       'optimcon() and grape_hilb() must handle a tiny Hilbert-space control problem.');

% Ensure that optimcon() has a process-pool ValueStore available
current_pool=gcp('nocreate');
if isempty(current_pool)
    parpool('Processes',1);
end

% Build a minimal Hilbert-space Spinach object for setup and GRAPE
spin_system=local_spin_system();

% Define one-spin operators and a two-step timing grid
S=pauli(2);
drift=sparse(2,2);
pulse_dt=[0.02 0.03];

% Configure a minimal optimal-control problem
control.operators={S.y};
control.rho_init={S.x};
control.rho_targ={S.z};
control.pwr_levels=1;
control.pulse_dt=pulse_dt;
control.drifts={{drift}};
control.method='lbfgs';
control.max_iter=0;
control.penalties={'none'};
control.p_weights=0;
control.l_bound=-100;
control.u_bound=100;
control.plotting={};
spin_system=optimcon(spin_system,control);

% Check that optimcon() absorbed the control metadata consistently
result=test_true(result,'optimcon control count',spin_system.control.ncontrols==1,...
                 'one supplied control operator must be registered');
result=test_close(result,'optimcon timing grid',spin_system.control.pulse_dt,pulse_dt,0,0,...
                  'the pulse timing grid must be stored without modification');
result=test_true(result,'optimcon rectangle grid',spin_system.control.pulse_ntpts==numel(pulse_dt),...
                 'rectangle integration must use one waveform value per interval');

% Evaluate the GRAPE fidelity and gradient for a non-trivial waveform
waveform=[7 11];
[traj_data,fidelity,grad]=grape_hilb(spin_system,{drift},control.operators,...
                                    waveform,S.x,S.z,'real');

% Build the independent exact Hilbert-space trajectory
rho=S.x;
for n=1:numel(pulse_dt)
    H=waveform(n)*S.y;
    P=expm(-1i*H*pulse_dt(n));
    rho=P*rho*P';
end
fid_ref=real(hdot(rho,S.z));

% Check the fidelity against direct matrix exponentiation
result=test_close(result,'grape_hilb fidelity',fidelity,fid_ref,1e-13,1e-13,...
                  'GRAPE fidelity must match direct Hilbert-space propagation');

% Compute a centred finite-difference gradient of the GRAPE fidelity
step_size=1e-6;
grad_ref=zeros(size(waveform));
for n=1:numel(waveform)
    wf_plus=waveform;
    wf_minus=waveform;
    wf_plus(n)=wf_plus(n)+step_size;
    wf_minus(n)=wf_minus(n)-step_size;
    [~,fid_plus]=grape_hilb(spin_system,{drift},control.operators,...
                            wf_plus,S.x,S.z,'real');
    [~,fid_minus]=grape_hilb(spin_system,{drift},control.operators,...
                             wf_minus,S.x,S.z,'real');
    grad_ref(n)=(fid_plus-fid_minus)/(2*step_size);
end

% Check the analytic GRAPE gradient against finite differences
result=test_close(result,'grape_hilb gradient',grad,grad_ref,1e-8,1e-8,...
                  'GRAPE gradient must match centred finite differences');

% Check that no heavy trajectory is returned when plotting is disabled
result=test_true(result,'grape_hilb trajectory suppression',isempty(traj_data.forward),...
                 'empty plotting requests should suppress trajectory storage');

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

