% Test: kernel/optimcon_weights | Weighted ensemble objectives
% Checks drift-weighted ensemble objectives and derivatives. Syntax:
%
%                  result=test_optimcon_weights()
%
% Outputs:
%
%     result - regression test result with explanatory messages
%
% Talos <talos@spindynamics.org>

function result=test_optimcon_weights()

% Build a quiet one-spin Hilbert-space problem
sys.isotopes={'1H'};
sys.magnet=1;
bas.formalism='zeeman-hilb';
bas.approximation={'none'};
spin_system=test_spin_system(sys,struct(),bas);
S=pauli(2);
waveform=[3 -1 2;2 4 -3];

% Correlate power with drift and replicate each drift over three offsets
control.isotopes={'1H'};
control.channels=[1;1];
control.operators={S.x,S.y};
control.rho_init={S.x};
control.rho_targ={S.y+0.3*S.z};
control.drifts={{-1.7*S.z},{2.9*S.z}};
control.pwr_levels=[0.8 1.2];
control.offsets={[-0.5 0.1 0.3]};
control.off_ops={S.z};
control.ens_corrs={'power_drift'};
control.phase_cycle=[0 0.7 0];
control.freeze=[true false false;false false false];
control.pulse_dt=[0.04 0.05 0.06];
control.method='newton';
control.max_iter=0;
control.plotting={};
control.traj_opts={'average'};
control.penalties={'NS'};
control.p_weights=0.2;
control.l_bound=-100;
control.u_bound=100;

% Evaluate each isolated drift, offset, and correlated power combination
single=cell(4,2);
raw_traj=cell(6,1);
for n=1:2

    % Fix the drift and its matched power
    local_control=control;
    local_control.drifts=control.drifts(n);
    local_control.pwr_levels=control.pwr_levels(n);
    local_control.traj_opts={};
    parts=cell(4,3);
    for k=1:3

        % Evaluate a single member, including the wrapper penalty channels
        local_control.offsets={control.offsets{1}(k)};
        local_system=optimcon(spin_system,local_control);
        local_system.control.return_traj=true;
        [traj,parts{2:4,k}]=grape_xy(waveform,local_system);
        raw_traj{n+2*(k-1)}=traj{1};
        parts{1,k}=traj{1}.forward;

    end

    % Average the three equally replicated offsets before applying drift weights
    for k=1:4
        single{k,n}=(1/3)*(parts{k,1}+parts{k,2}+parts{k,3});
    end

end

% Test an unequal quadrature, a uniform quadrature, and a zero-weight member
result=new_test_result('kernel/optimcon_weights','Weighted ensemble objectives',...
                       'Checks case-order weighting after physical waveform transformations.');
result=test_true(result,'distinct objectives',abs(single{2,1}(1)-single{2,2}(1))>1e-3,...
                 'The individual drift fidelities must distinguish the two weights.');
for weights={[1 3],[1 1],[1 0]}

    % Construct the ensemble and compare its physical and penalty channels
    local_control=control;
    local_control.drift_weights=weights{1};
    local_system=optimcon(spin_system,local_control);
    local_system.control.return_traj=true;
    [traj,fid,grad,hess]=grape_xy(waveform,local_system);
    w=weights{1}/sum(weights{1});
    label=sprintf('%g:%g',weights{1});
    result=test_close(result,[label ' fidelity'],fid,w(1)*single{2,1}+w(2)*single{2,2},1e-12,0,...
                      'Weighted physical fidelities and unchanged penalty channels.');
    result=test_close(result,[label ' gradient'],grad,w(1)*single{3,1}+w(2)*single{3,2},1e-12,0,...
                      'The weights multiply the pulled-back gradients.');
    result=test_close(result,[label ' Hessian'],hess,w(1)*single{4,1}+w(2)*single{4,2},1e-12,0,...
                      'The weights multiply the transformed Hessians.');
    expected=w(1)*single{1,1}+w(2)*single{1,2};
    result=test_close(result,[label ' trajectory'],horzcat(traj{1}.forward{:}),...
                      horzcat(expected{:}),1e-12,0,'Average trajectory weights match the objective.');

end

% Verify reordered diagnostic trajectories without applying case weights
local_control.drift_weights=[1 3];
local_control.traj_opts={};
local_system=optimcon(spin_system,local_control);
local_system.control.return_traj=true;
order=[local_system.control.worker_cases{:}];
result=test_true(result,'nontrivial worker ordering',~isequal(order,1:6),...
                 'The fixture must exercise drift sorting on a non-drift-major catalog.');
[traj,~]=grape_xy(waveform,local_system);
for n=1:numel(traj)

    % Check the unweighted trajectory in catalog order
    result=test_close(result,sprintf('raw trajectory %d',n),horzcat(traj{n}.forward{:}),...
                      horzcat(raw_traj{n}.forward{:}),1e-12,0,'Raw trajectories retain catalog order and scale.');

end

% Check that explicit uniform weights reproduce the default
uniform=control;
uniform.drift_weights=[1 1];
[~,fid_u,grad_u,hess_u]=grape_xy(waveform,optimcon(spin_system,uniform));
[~,fid_d,grad_d,hess_d]=grape_xy(waveform,optimcon(spin_system,control));
result=test_close(result,'default fidelity',fid_u,fid_d,0,0,'The default gives equal weight to all retained cases.');
result=test_close(result,'default gradient',grad_u,grad_d,0,0,'Uniform explicit and default gradients agree.');
result=test_close(result,'default Hessian',hess_u,hess_d,0,0,'Uniform explicit and default Hessians agree.');

% Refuse stochastic budgets and all-zero weights
for n=1:2

    % Select the invalid input and expected error
    invalid=control;
    if n==1
        invalid.drift_weights=[1 3];
        invalid.budget=0.5;
        err_id='Spinach:optimcon:weightedBudget';
    else
        invalid.drift_weights=[0 0];
        err_id='Spinach:optimcon:zeroDriftWeight';
    end

    % Check the public setup rejection
    refused=false;
    try
        optimcon(spin_system,invalid);
    catch err
        refused=strcmp(err.identifier,err_id);
    end
    result=test_true(result,err_id,refused,'Invalid weighted quadrature must fail at setup.');

end

end


