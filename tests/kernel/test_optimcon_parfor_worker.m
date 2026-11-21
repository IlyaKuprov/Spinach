% Tests optimal control launched inside an external parfor loop. Syntax:
%
%                   result=test_optimcon_parfor_worker()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test optimises several distinct one-spin Hilbert-space problems,
% some with power-level ensembles, with optimcon() and fmaxnewton()
% called inside a parfor body, and compares each worker waveform with
% the same optimisation on the client. It also checks that ensemble()
% rejects a problem prepared by optimcon() on the client pool when it
% is evaluated on a worker.
%
% ilya.kuprov@weizmann.ac.il

function result=test_optimcon_parfor_worker()

% Announce the test target
fprintf('TESTING: Optimal control inside an external parfor loop\n');

% State the optimal-control target of the test
result=new_test_result('optimcon/parfor_worker',...
                       'Optimal control inside an external parfor loop',...
                       'optimcon() and fmaxnewton() must run inside parfor and reproduce client optimisations.');

% Ensure that a parallel pool exists for the external loop
current_pool=gcp('nocreate');
if isempty(current_pool)
    parpool('Processes',2);
end

% Distinct problems: drift offsets and power-level ensembles
offsets={0 0.5 1.0 1.5}; pwr_levels={1 [0.9 1.0 1.1] 1 [0.95 1.05]};
n_probs=numel(offsets);

% Optimise every problem on the client
x_client=cell(1,n_probs);
for n=1:n_probs
    x_client{n}=local_optimise(offsets{n},pwr_levels{n});
end

% Optimise the same problems inside parfor
x_worker=cell(1,n_probs); messages=cell(1,n_probs);
parfor n=1:n_probs
    try
        x_worker{n}=local_optimise(offsets{n},pwr_levels{n}); messages{n}='';
    catch err
        messages{n}=err.message;
    end
end

% Compare worker and client waveforms
for n=1:n_probs
    label=['parfor problem ' int2str(n)];
    result=test_true(result,[label ' optimisation'],isempty(messages{n}),...
                     strtrim(['the optimisation must complete on a worker. ' messages{n}]));
    if isempty(messages{n})
        result=test_close(result,[label ' waveform'],x_worker{n},x_client{n},0,1e-10,...
                          'the worker optimisation must reproduce the client optimisation.');
    end
end

% Prepare a problem on the client pool
[spin_system,control,guess]=local_problem(0,1);
spin_system=optimcon(spin_system,control);

% Evaluate the client-prepared problem on workers
rejected=false(1,2);
parfor n=1:2
    try
        [~,~]=grape_xy(guess,spin_system);
    catch err
        rejected(n)=contains(err.message,'optimcon()');
    end
end
result=test_true(result,'client-prepared problem on a worker',all(rejected),...
                 'ensemble() must reject a problem prepared on a different pool.');

end

% Optimises one test problem
function x=local_optimise(offset,pwr_levels)

% Set up and run the optimisation
[spin_system,control,guess]=local_problem(offset,pwr_levels);
spin_system=optimcon(spin_system,control);
x=fmaxnewton(spin_system,@grape_xy,guess);

end

% Builds a one-spin state transfer problem with a frequency offset
function [spin_system,control,guess]=local_problem(offset,pwr_levels)

% Get the spin system and Pauli matrices
spin_system=local_spin_system(); S=pauli(2);

% Transfer Z magnetisation into X magnetisation under an offset
control.isotopes={'E'}; control.channels=[1; 1];
control.drifts={{offset*S.z}}; control.operators={S.x,S.y};
control.rho_init={S.z/norm(full(S.z),'fro')};
control.rho_targ={S.x/norm(full(S.x),'fro')};

% Short L-BFGS optimisation on a fixed grid
control.pwr_levels=pwr_levels; control.pulse_dt=0.1*ones(1,20);
control.method='lbfgs'; control.max_iter=10;
control.penalties={'none'}; control.p_weights=0;
control.l_bound=-100; control.u_bound=100; control.plotting={};

% Deterministic initial guess
t=1:20; guess=0.1*[cos(t); sin(t)];

end

% Builds a minimal quiet Spinach object for Hilbert-space helpers
function spin_system=local_spin_system()

% Set quiet output and tight tolerances
spin_system.sys.output='hush';
spin_system.sys.enable={};
spin_system.sys.disable={};
spin_system.tols.liouv_zero=1e-14;
spin_system.tols.small_matrix=64;
spin_system.tols.dense_matrix=0.5;
spin_system.tols.prop_chop=1e-14;
spin_system.bas.formalism='zeeman-hilb';
spin_system.comp.isotopes={'E'};

end

