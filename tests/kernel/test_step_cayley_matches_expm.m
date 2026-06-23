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

end

