% Tests fourth-order convergence of the Cayley-Magnus stepper. Syntax:
%
%                    result=test_step_cayley_order()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% ilya.kuprov@weizmann.ac.il

function result=test_step_cayley_order()

% Announce the test target
fprintf('TESTING: Cayley-Magnus fourth-order convergence\n');

% State the propagation target of the test
result=new_test_result('kernel/step_cayley_order',...
                       'Cayley-Magnus fourth-order convergence',...
                       'step_cayley() must show fourth-order convergence.');

% Build a one-proton Liouville-space spin system
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='zeeman-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Define a non-commuting constant generator
S=pauli(2); L=2*pi*(S.x+0.7*S.y+0.2*S.z);
rho=[1; -0.5]; dt=0.05;

% Compare two Cayley-Magnus refinements against expm
rho_ref=expm(-1i*L*dt)*rho;
err_1=norm(step_cayley(spin_system,L,rho,dt)-rho_ref,2);

% Take two half steps to reach the same final time
rho_half=step_cayley(spin_system,L,rho,dt/2);
rho_half=step_cayley(spin_system,L,rho_half,dt/2);
err_2=norm(rho_half-rho_ref,2);

% Check that halving the step improves by fourth-order margin
result=test_true(result,'fourth-order error reduction',err_1/err_2>10,...
                 'halving the time step should strongly reduce Cayley-Magnus error');

end

