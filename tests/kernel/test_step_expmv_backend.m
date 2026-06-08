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

end
