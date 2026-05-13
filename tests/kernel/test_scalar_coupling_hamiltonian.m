% Tests the two-spin scalar-coupling Hamiltonian. Syntax:
%
%                    result=test_scalar_coupling_hamiltonian()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that an isotropic scalar coupling J produces the textbook
% Hamiltonian 2*pi*J*(Ix*Sx+Iy*Sy+Iz*Sz).
%
% ilya.kuprov@weizmann.ac.il

function result=test_scalar_coupling_hamiltonian()

% Announce the test target
fprintf('TESTING: Scalar coupling Hamiltonian\n');

% State the Hamiltonian target of the test
result=new_test_result('kernel/scalar_coupling_hamiltonian',...
                       'Scalar coupling Hamiltonian',...
                       'an isotropic J coupling must produce 2*pi*J I dot S.');

% Build a two-proton Hilbert-space spin system with a 10 Hz J coupling
sys.magnet=0;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0,0};
inter.coupling.scalar{1,2}=10;
inter.coupling.scalar{2,2}=0;
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Build Spinach and textbook Hamiltonians
H_obs=hamiltonian(assume(spin_system,'nmr'));
IxSx=operator(spin_system,{'Lx','Lx'},{1,2});
IySy=operator(spin_system,{'Ly','Ly'},{1,2});
IzSz=operator(spin_system,{'Lz','Lz'},{1,2});
H_ref=2*pi*10*(IxSx+IySy+IzSz);

% Check the scalar-coupling Hamiltonian
result=test_close(result,'isotropic J Hamiltonian',H_obs,H_ref,1e-9,1e-12,...
                  'scalar coupling is rotationally invariant I dot S in rad/s units');

end

