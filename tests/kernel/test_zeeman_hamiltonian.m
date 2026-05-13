% Tests the one-spin Zeeman Hamiltonian. Syntax:
%
%                    result=test_zeeman_hamiltonian()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks Spinach's NMR convention for a positive chemical shift:
% the rotating-frame Hamiltonian contribution is -2*pi*nu*Lz.
%
% ilya.kuprov@weizmann.ac.il

function result=test_zeeman_hamiltonian()

% Announce the test target
fprintf('TESTING: Zeeman Hamiltonian sign and units\n');

% State the Hamiltonian target of the test
result=new_test_result('kernel/zeeman_hamiltonian',...
                       'Zeeman Hamiltonian sign and units',...
                       'a scalar chemical shift must enter the NMR Hamiltonian with Spinach sign and rad/s units.');

% Build a one-proton Hilbert-space spin system with a 1 ppm shift
sys.magnet=14.1;
sys.isotopes={'1H'};
inter.zeeman.scalar={1};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Build Spinach and reference Hamiltonians
H_obs=hamiltonian(assume(spin_system,'nmr'));
nu=ppm2hz(1,sys.magnet,'1H');
H_ref=-2*pi*nu*operator(spin_system,'Lz',1);

% Check the physical frequency and sign convention
result=test_close(result,'one-spin Zeeman Hamiltonian',H_obs,H_ref,1e-6,1e-12,...
                  'positive ppm gives -omega*Lz in the Spinach NMR rotating-frame convention');

end

