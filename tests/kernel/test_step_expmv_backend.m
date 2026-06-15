% Tests step() exponential-action backend selection. Syntax:
%
%                    result=test_step_expmv_backend()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that sys.expmv_backend='auto' matches the native
% reordered Taylor step path and that unsupported step() modes fall back
% to the default backend.
%
% ilya.kuprov@weizmann.ac.il

function result=test_step_expmv_backend()

% Announce the test target
fprintf('TESTING: step() exponential-action backend selection\n');

% State the propagation target of the test
result=new_test_result('kernel/step_expmv_backend',...
                       'step() exponential-action backend selection',...
                       'automatic step() backend selection must match the native Taylor path.');

% Build a one-spin Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={3.0};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);
spin_system.tols.small_matrix=2;

% Define a simple state and generator
L=operator(spin_system,'Lz','1H');
rho=state(spin_system,'Lx','1H')+0.25*state(spin_system,'Ly','1H');
dt=2.5e-4;

% Build automatic and default systems
spin_auto=spin_system;
spin_auto.sys.expmv_backend='auto';
spin_default=spin_system;
spin_default.sys.expmv_backend='default';

% Vector stacks are intentionally routed through the default stack path
rho_stack=[rho 2*rho];
rho_obs=step(spin_auto,L,rho_stack,dt);
rho_ref=step(spin_default,L,rho_stack,dt);
result=test_close(result,'step stack backend fallback',rho_obs,rho_ref,1e-12,1e-12,...
                  'automatic backend selection must preserve stack propagation');

% High-norm single-vector propagation exercises the automatic selector
fast_dt=10/cheap_norm(L);
rho_obs=step(spin_auto,L,rho,fast_dt);
rho_ref=step(spin_default,L,rho,fast_dt);
result=test_close(result,'step auto backend single vector',rho_obs,rho_ref,1e-12,1e-12,...
                  'automatic single-vector backend selection must match the native Taylor path');

% Dense CPU problems in this range are routed through Matlab expmv()
dense_dim=128;
dense_op=diag(linspace(-6,6,dense_dim));
dense_rho=sin((1:dense_dim).'/dense_dim);
stats=struct('matrix',dense_op,'time_step',1,'dimension',dense_dim,...
             'is_sparse',false,'is_gpu',false,'norm_mat',6);
backends=struct('default',@(~)'default','expmv',@(~)'expmv',...
                'tay1',@(~)'tay1','tay2',@(~)'tay2');
backend=step_heuristics(stats,backends);
result=test_true(result,'step dense expmv selection',strcmp(backend([]),'expmv'),...
                 'dense CPU heuristic must select Matlab expmv in its documented range');
rho_obs=step(spin_auto,dense_op,dense_rho,1);
rho_ref=step(spin_default,dense_op,dense_rho,1);
result=test_close(result,'step dense expmv backend',rho_obs,rho_ref,1e-11,1e-11,...
                  'Matlab expmv backend must match the native Taylor path');

% GPU Taylor backends must keep scalar policy values on the CPU
if (exist('gpuDeviceCount','file')==2)&&(gpuDeviceCount>0)
    spin_gpu=spin_auto; spin_gpu.sys.enable={'gpu'};
    spin_gpu_def=spin_default; spin_gpu_def.sys.enable={'gpu'};
    gpu_dim=64;
    gpu_rho=cos((1:gpu_dim).'/gpu_dim);
    for alpha=[6 20]
        gpu_op=diag(linspace(-alpha,alpha,gpu_dim));
        rho_obs=gather(step(spin_gpu,gpu_op,gpu_rho,1));
        rho_ref=gather(step(spin_gpu_def,gpu_op,gpu_rho,1));
        result=test_close(result,['step GPU auto backend alpha=' num2str(alpha)],...
                          rho_obs,rho_ref,1e-10,1e-10,...
                          'GPU automatic Taylor backends must match the native Taylor path');
    end
end

% Product quadrature inputs are not eligible for the automatic selector
Lx=operator(spin_system,'Lx','1H');
rho_obs=step(spin_auto,{L+0.5*Lx,L-0.25*Lx},rho,dt);
rho_ref=step(spin_default,{L+0.5*Lx,L-0.25*Lx},rho,dt);
result=test_close(result,'step two-point backend fallback',rho_obs,rho_ref,1e-12,1e-12,...
                  'automatic backend selection must preserve two-point quadrature');

rho_obs=step(spin_auto,{L+0.5*Lx,L+0.1*Lx,L-0.25*Lx},rho,dt);
rho_ref=step(spin_default,{L+0.5*Lx,L+0.1*Lx,L-0.25*Lx},rho,dt);
result=test_close(result,'step three-point backend fallback',rho_obs,rho_ref,1e-12,1e-12,...
                  'automatic backend selection must preserve three-point quadrature');

% Complex time is deliberately excluded from automatic backend selection
rho_obs=step(spin_auto,L,rho,1i*dt);
rho_ref=step(spin_default,L,rho,1i*dt);
result=test_close(result,'step complex-time backend fallback',rho_obs,rho_ref,1e-12,1e-12,...
                  'automatic backend selection must preserve complex-time propagation');

% Malformed heuristic inputs must trigger Spinach validation
heuristic_err=false();
try
    step_heuristics(struct(),struct());
catch err
    heuristic_err=contains(err.message,'missing required fields');
end

result=test_true(result,'step heuristics grumble',heuristic_err,...
                 'malformed heuristic inputs must produce a clear validation error');

% Large-alpha scalar propagation must remain stable
spin_auto.sys.output='hush';
spin_default.sys.output='hush';
rho_obs=step(spin_auto,1000,1,1);
rho_ref=step(spin_default,1000,1,1);
result=test_close(result,'step large-alpha scalar backend',rho_obs,rho_ref,1e-10,1e-10,...
                  'automatic backend selection must remain stable for large scalar generators');

% Sparse automatic Taylor backend must avoid dense history arrays
sparse_dim=20000;
sparse_op=spdiags([ones(sparse_dim,1) -2*ones(sparse_dim,1) ones(sparse_dim,1)],...
                 -1:1,sparse_dim,sparse_dim);
sparse_op=(40/cheap_norm(sparse_op))*sparse_op;
sparse_rho=ones(sparse_dim,1);
profile clear
profile('-memory','on');
step(spin_auto,sparse_op,sparse_rho,1);
profile off
prof_info=profile('info');
prof_names={prof_info.FunctionTable.FunctionName};
tay_rows=contains(prof_names,'step_expmv_tay');
if any(tay_rows)
    sparse_peak=max([prof_info.FunctionTable(tay_rows).PeakMem]);
else
    sparse_peak=Inf;
end
dense_hist=16*sparse_dim*61;
result=test_true(result,'step sparse rolling storage',sparse_peak<dense_hist/2,...
                 'sparse automatic Taylor backends must not allocate dense history arrays');

end
