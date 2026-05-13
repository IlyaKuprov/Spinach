% Tests dynamic propagation front-end kernels on tiny systems. Syntax:
%
%                    result=test_dynamic_propagation_frontends()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test exercises propagator(), step(), evolution(), krylov(), and
% reduce() against direct finite-dimensional propagation references.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_propagation_frontends()

% Announce the test target
fprintf('TESTING: Dynamic propagation front ends\n');

% State the dynamic propagation target of the test
result=new_test_result('kernel/dynamic_propagation_frontends',...
                       'Dynamic propagation front ends',...
                       'propagation front ends must match direct finite-dimensional references.');

% Check scaled Taylor propagator and step() branches
result=local_test_propagator_step(result);

% Check evolution() output modes against explicit propagator products
result=local_test_evolution(result);

% Check direct krylov() output modes against step() references
result=local_test_krylov(result);

% Check reduce() projector invariants and blanket-disable branch
result=local_test_reduce(result);

end


function result=local_test_propagator_step(result)

% Build a one-spin Liouville-space system and force Taylor propagation
spin_system=local_liouville_system();
spin_system.tols.small_matrix=2;
L=operator(spin_system,'Lz','1H');
rho=state(spin_system,'Lx','1H')+0.25*state(spin_system,'Ly','1H');
dt=2.5e-4;

% Compare propagator() against direct matrix exponentiation
P_obs=propagator(spin_system,L,dt);
P_ref=expm(full(-1i*L*dt));
result=test_close(result,'propagator Taylor branch',P_obs,P_ref,1e-10,1e-10,...
                  'propagator() must reproduce expm() to the configured sparse chop tolerance');

% Compare numeric step() against the same propagator action
rho_ref=P_ref*rho;
rho_obs=step(spin_system,L,rho,dt);
result=test_close(result,'step numeric Liouville branch',rho_obs,rho_ref,1e-12,1e-12,...
                  'step() must propagate Liouville-space vectors with the Spinach sign convention');

% Check the zero-time shortcut
rho_obs=step(spin_system,L,rho,0);
result=test_close(result,'step zero-time branch',rho_obs,rho,1e-14,1e-14,...
                  'step() must return the input state unchanged for a zero time step');

% Check two-point and three-point quadrature paths in the constant limit
rho_obs=step(spin_system,{L,L},rho,dt);
result=test_close(result,'step two-point constant branch',rho_obs,rho_ref,1e-12,1e-12,...
                  'step() two-point quadrature must reduce to the constant-generator result');
rho_obs=step(spin_system,{L,L,L},rho,dt);
result=test_close(result,'step three-point constant branch',rho_obs,rho_ref,1e-12,1e-12,...
                  'step() three-point quadrature must reduce to the constant-generator result');

% Check the state-dependent generator dispatch with a constant callback
rho_obs=step(spin_system,{@(time,state)L,0.0,'PWCL'},rho,dt);
result=test_close(result,'step function-handle branch',rho_obs,rho_ref,1e-12,1e-12,...
                  'step() must route state-dependent generator callbacks through iserstep()');

% Check the Hilbert-space commutator-series branch against unitary evolution
spin_system=local_hilbert_system();
spin_system.tols.small_matrix=1;
S=pauli(2);
H=2*pi*123*S.z;
rho_h=S.x+0.25*S.y;
P=expm(-1i*H*dt);
rho_ref=P*rho_h*P';
rho_obs=step(spin_system,H,rho_h,dt);
result=test_close(result,'step Hilbert commutator branch',rho_obs,rho_ref,1e-12,1e-12,...
                  'step() must propagate Hilbert-space density matrices by unitary similarity');

end


function result=local_test_evolution(result)

% Build a one-spin Liouville-space system with trajectory reduction disabled
spin_system=local_liouville_system();
spin_system.sys.disable=unique([spin_system.sys.disable {'trajlevel'}]);
L=operator(spin_system,'Lz','1H');
rho=state(spin_system,'Lx','1H')+0.25*state(spin_system,'Ly','1H');
coil_x=state(spin_system,'Lx','1H');
coil_y=state(spin_system,'Ly','1H');
dt=2.5e-4;
nsteps=3;
P=propagator(spin_system,L,dt);
traj_ref=[rho P*rho P^2*rho P^3*rho];

% Check final-state evolution
rho_obs=evolution(spin_system,L,[],rho,dt,nsteps,'final');
result=test_close(result,'evolution final output',rho_obs,traj_ref(:,end),1e-12,1e-12,...
                  'evolution() final mode must match repeated propagator action');

% Check trajectory evolution
traj_obs=evolution(spin_system,L,[],rho,dt,nsteps,'trajectory');
result=test_close(result,'evolution trajectory output',traj_obs,traj_ref,1e-12,1e-12,...
                  'evolution() trajectory mode must return every propagated state vector');

% Check single-channel observable evolution
fid_obs=evolution(spin_system,L,coil_x,rho,dt,nsteps,'observable');
fid_ref=coil_x'*traj_ref;
result=test_close(result,'evolution observable output',fid_obs,fid_ref.',1e-12,1e-12,...
                  'evolution() observable mode must detect every trajectory point');

% Check multichannel observable evolution
fid_obs=evolution(spin_system,L,[coil_x coil_y],rho,dt,nsteps,'multichannel');
fid_ref=[coil_x';coil_y']*traj_ref;
result=test_close(result,'evolution multichannel output',fid_obs,fid_ref,1e-12,1e-12,...
                  'evolution() multichannel mode must detect every coil row separately');

% Check refocusing evolution over a stack of initial states
rho_stack=[rho 2*rho 3*rho];
rho_ref=[rho_stack(:,1) P*rho_stack(:,2) P^2*rho_stack(:,3)];
rho_obs=evolution(spin_system,L,[],rho_stack,dt,2,'refocus');
result=test_close(result,'evolution refocus output',rho_obs,rho_ref,1e-12,1e-12,...
                  'evolution() refocus mode must propagate the nth stack member for n-1 steps');

% Check total integral mode with a scalar damping generator
L_decay=-1i*speye(size(L,1));
total_obs=evolution(spin_system,L_decay,coil_x,rho,dt,nsteps,'total');
total_ref=real(coil_x'*rho);
result=test_close(result,'evolution total output',total_obs,total_ref,1e-12,1e-12,...
                  'evolution() total mode must compute the infinite-time damped observable integral');

end


function result=local_test_krylov(result)

% Build a one-spin Liouville-space system for direct Krylov propagation
spin_system=local_liouville_system();
L=operator(spin_system,'Lz','1H');
rho=state(spin_system,'Lx','1H')+0.25*state(spin_system,'Ly','1H');
coil_x=state(spin_system,'Lx','1H');
coil_y=state(spin_system,'Ly','1H');
dt=2.5e-4;
nsteps=3;
P=propagator(spin_system,L,dt);
traj_ref=[rho P*rho P^2*rho P^3*rho];

% Check final and trajectory Krylov modes
rho_obs=krylov(spin_system,L,[],rho,dt,nsteps,'final');
result=test_close(result,'krylov final output',rho_obs,traj_ref(:,end),1e-12,1e-12,...
                  'krylov() final mode must match repeated step() propagation');
traj_obs=krylov(spin_system,L,[],rho,dt,nsteps,'trajectory');
result=test_close(result,'krylov trajectory output',traj_obs,traj_ref,1e-12,1e-12,...
                  'krylov() trajectory mode must return every propagated state vector');

% Check observable and multichannel Krylov modes
fid_obs=krylov(spin_system,L,coil_x,rho,dt,nsteps,'observable');
fid_ref=coil_x'*traj_ref;
result=test_close(result,'krylov observable output',fid_obs,fid_ref.',1e-12,1e-12,...
                  'krylov() observable mode must detect every trajectory point');
fid_obs=krylov(spin_system,L,[coil_x coil_y],rho,dt,nsteps,'multichannel');
fid_ref=[coil_x';coil_y']*traj_ref;
result=test_close(result,'krylov multichannel output',fid_obs,fid_ref,1e-12,1e-12,...
                  'krylov() multichannel mode must detect every coil row separately');

% Check refocus Krylov mode over a stack of initial states
rho_stack=[rho 2*rho 3*rho];
rho_ref=[rho_stack(:,1) P*rho_stack(:,2) P^2*rho_stack(:,3)];
rho_obs=krylov(spin_system,L,[],rho_stack,dt,[], 'refocus');
result=test_close(result,'krylov refocus output',rho_obs,rho_ref,1e-12,1e-12,...
                  'krylov() refocus mode must propagate the nth stack member for n-1 steps');

end


function result=local_test_reduce(result)

% Build a one-spin Liouville-space system for projector checks
spin_system=local_liouville_system();
L=operator(spin_system,'Lz','1H');
rho=state(spin_system,'Lx','1H');
projectors=reduce(spin_system,L,rho);

% Check that every returned projector has orthonormal columns
for n=1:numel(projectors)
    P=projectors{n};
    result=test_close(result,['reduce projector orthogonality ' num2str(n)],P'*P,speye(size(P,2)),1e-14,1e-14,...
                      'reduce() projectors must be column-orthonormal');
end

% Check that the active initial state is preserved by all projectors together
rho_rec=sparse(size(rho,1),1);
for n=1:numel(projectors)
    P=projectors{n};
    rho_rec=rho_rec+P*(P'*rho);
end
result=test_close(result,'reduce state preservation',rho_rec,rho,1e-14,1e-14,...
                  'reduce() must retain the source state in the projected subspaces');

% Check the blanket trajectory-level disable branch
spin_system.sys.disable=unique([spin_system.sys.disable {'trajlevel'}]);
projectors=reduce(spin_system,L,rho);
result=test_true(result,'reduce trajlevel disable branch',isscalar(projectors)&&isequal(projectors{1},1),...
                 'reduce() must return the scalar identity when trajectory-level reduction is disabled');

end


function spin_system=local_liouville_system()

% Build the common one-spin spherical-tensor Liouville-space system
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={0.0};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);
spin_system=assume(spin_system,'nmr');

end


function spin_system=local_hilbert_system()

% Build the common one-spin Hilbert-space system
sys.magnet=0.0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0.0};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

end


