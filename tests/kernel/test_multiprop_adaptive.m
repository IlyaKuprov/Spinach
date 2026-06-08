% Tests adaptive repeated propagator application. Syntax:
%
%                    result=test_multiprop_adaptive()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks binary adaptive squaring in multiprop() against explicit
% matrix-power references for state vectors and density matrices.
%
% ilya.kuprov@weizmann.ac.il

function result=test_multiprop_adaptive()

% Announce the test target
fprintf('TESTING: Adaptive repeated propagator application\n');

% State the propagation target of the test
result=new_test_result('kernel/multiprop_adaptive',...
                       'Adaptive repeated propagator application',...
                       'multiprop() must apply propagators repeatedly by binary adaptive squaring.');

% Define a non-normal sparse propagator and a state vector
P=sparse([0.90 0.10 0.00;0.00 0.80 0.20;0.05 0.00 0.70]);
rho=[1.0;2.0;3.0];
N=uint64(13);

% Compare state-vector propagation with an explicit matrix power
rho_obs=multiprop(P,rho,N);
rho_ref=(P^double(N))*rho;
result=test_close(result,'state vector propagation',rho_obs,rho_ref,1e-14,1e-14,...
                  'Liouville-space and wavefunction states are propagated by left multiplication');

% Check the zero-step shortcut
rho_obs=multiprop(P,rho,0);
result=test_close(result,'zero-step propagation',rho_obs,rho,1e-15,1e-15,...
                  'zero applications of a propagator leave the state unchanged');

% Check scalar vector branch consistency
rho_obs=multiprop(0.5,2.0,3);
result=test_close(result,'scalar vector propagation',rho_obs,0.25,1e-15,1e-15,...
                  'one-dimensional state vectors follow the vector propagation branch');

% Define a Hilbert-space unitary propagator and a density matrix
S=pauli(2);
H=2*pi*(0.30*S.x+0.70*S.z);
P=expm(-1i*H*0.125);
rho=S.x+0.25*S.y+0.10*S.z;
N=9;

% Compare density-matrix propagation with an explicit matrix power
P_ref=P^N;
rho_obs=multiprop(P,rho,N);
rho_ref=P_ref*rho*P_ref';
result=test_close(result,'density matrix propagation',rho_obs,rho_ref,1e-13,1e-13,...
                  'Hilbert-space density matrices are propagated by left and right multiplication');

% Check row-vector rejection
caught=false;
try
    multiprop(P,[1 0],3);
catch
    caught=true;
end

result=test_true(result,'row vector rejection',caught,...
                 'row vectors are not valid Spinach state vectors');

end

