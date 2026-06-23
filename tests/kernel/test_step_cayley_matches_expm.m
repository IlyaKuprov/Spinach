% Tests Cayley-Magnus vector propagation against matrix exponentiation.
% Syntax:
%
%                    result=test_step_cayley_matches_expm()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% ilya.kuprov@weizmann.ac.il

function result=test_step_cayley_matches_expm()

% Announce the test target
fprintf('TESTING: Cayley-Magnus vector propagation against expm\n');

% State the propagation target of the test
result=new_test_result('kernel/step_cayley_matches_expm',...
                       'Cayley-Magnus vector propagation against expm',...
                       'step_cayley() must reproduce constant-generator vector propagation.');

% Build a one-proton Liouville-space spin system
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='zeeman-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Define a generator and an initial state stack
S=pauli(2); L=kron(S.z,eye(2))-kron(eye(2),transpose(S.z));
rho=[1; 0.25; -0.5; 0.75];
rho=[rho 0.5*rho+0.1i*flipud(rho)];
dt=2.5e-3;

% Build the independent exact propagator
rho_ref=expm(-1i*L*dt)*rho;
rho_obs=step_cayley(spin_system,L,rho,dt);

% Check vector propagation
result=test_close(result,'step_cayley versus expm',rho_obs,rho_ref,1e-10,1e-10,...
                  'constant-generator vector propagation matches matrix exponentiation');

% Increase the step to exercise norm-estimated Cayley subdivision
dt=0.25;

% Build the independent exact propagator
rho_ref=expm(-1i*L*dt)*rho;
rho_obs=step_cayley(spin_system,L,rho,dt);

% Check scaled vector propagation
result=test_close(result,'step_cayley norm-estimated scaling',rho_obs,rho_ref,1e-6,1e-6,...
                  'norm-estimated Cayley subdivision controls the rational defect');

% Check non-finite time step rejection
try
    step_cayley(spin_system,L,rho,Inf);
    finite_rejected=false();
catch err
    finite_rejected=contains(err.message,'finite scalar');
end
result=test_true(result,'non-finite time step rejection',finite_rejected,...
                 'non-finite time steps must be rejected at input validation');

% Define a linearly varying non-commuting generator
S=pauli(2);
dt=0.2; nsteps=400;
L_left=12*(S.x+0.7*S.y);
L_right=12*(S.x+0.7*S.y+0.9*S.z);
rho=[1; -0.3+0.2i];

% Build the midpoint product reference
rho_ref=rho; time_slice=dt/nsteps;
for n=1:nsteps

    % Evaluate the linearly interpolated generator
    time_pos=(n-1/2)*time_slice;
    L_inst=L_left+(time_pos/dt)*(L_right-L_left);

    % Apply the local midpoint exponential
    rho_ref=expm(-1i*L_inst*time_slice)*rho_ref;

end

% Propagate from endpoint samples
rho_obs=step_cayley(spin_system,{L_left,L_right},rho,dt);

% Check sampled time-dependent propagation
result=test_close(result,'step_cayley endpoint samples',rho_obs,rho_ref,3e-2,3e-2,...
                  'endpoint-sampled time-dependent propagation matches midpoint reference');

end
