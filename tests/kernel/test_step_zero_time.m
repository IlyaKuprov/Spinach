% Tests zero-duration propagation. Syntax:
%
%                    result=test_step_zero_time()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks the identity limit of the propagator: a zero time step
% must leave the density matrix exactly unchanged.
%
% ilya.kuprov@weizmann.ac.il

function result=test_step_zero_time()

% State the propagation target of the test
result=new_test_result('kernel/step_zero_time',...
                       'Zero-duration propagation identity',...
                       'a propagator over zero time is the identity map.');

% Build a one-proton Hilbert-space spin system
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Propagate an arbitrary Hermitian density matrix for zero time
S=pauli(2);
rho=S.x+2*S.z;
H=3*S.x+5*S.z;
rho_obs=step(spin_system,H,rho,0);

% Check the identity limit
result=test_close(result,'rho(t=0)=rho(0)',rho_obs,rho,1e-15,1e-15,...
                  'zero-duration evolution cannot change the state');

end

