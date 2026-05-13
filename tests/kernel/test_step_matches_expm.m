% Tests Hilbert-space propagation against matrix exponentiation. Syntax:
%
%                    result=test_step_matches_expm()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks the Spinach sign convention for density-matrix evolution:
% rho(t)=exp(-iHt) rho(0) exp(+iHt).
%
% ilya.kuprov@weizmann.ac.il

function result=test_step_matches_expm()

% Announce the test target
fprintf('TESTING: Hilbert propagation against expm\n');

% State the propagation target of the test
result=new_test_result('kernel/step_matches_expm',...
                       'Hilbert propagation against expm',...
                       'step() must reproduce unitary density-matrix propagation.');

% Build a one-proton Hilbert-space spin system
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Define a Hamiltonian and an initial density matrix
S=pauli(2);
H=2*pi*123*S.z;
rho=S.x+0.25*S.y;
dt=2.5e-3;

% Build the independent exact propagator
P=expm(-1i*H*dt);
rho_ref=P*rho*P';
rho_obs=step(spin_system,H,rho,dt);

% Check exact finite-dimensional propagation
result=test_close(result,'step versus expm',rho_obs,rho_ref,1e-13,1e-13,...
                  'finite Hilbert-space propagation is exactly unitary');

end

