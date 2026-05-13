% Tests thermal equilibrium state construction paths. Syntax:
%
%              result=test_dynamic_state_equilibrium_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks Hilbert-space, Zeeman-Liouville, and oriented-Hamiltonian
% thermal equilibrium construction against direct Boltzmann references.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_state_equilibrium_suite()

% Announce the test target
fprintf('TESTING: Thermal equilibrium constructors\n');

% State the equilibrium-constructor target of the test
result=new_test_result('kernel/dynamic_state_equilibrium_suite',...
                       'Thermal equilibrium constructors',...
                       'equilibrium() must produce normalised Boltzmann states in supported formalisms.');

% Build a one-spin Hilbert-space system with finite temperature
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
inter.temperature=300;
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_h=test_spin_system(sys,inter,bas);

% Set an explicit non-degenerate Hamiltonian in angular frequency units
H=2*pi*diag([-5e11 3e11]);
beta=spin_h.tols.hbar/(spin_h.tols.kbol*spin_h.rlx.temperature);
rho_ref=expm(-beta*H);
rho_ref=rho_ref/trace(rho_ref);

% Compare Hilbert-space equilibrium to the direct Boltzmann density matrix
rho_h=equilibrium(spin_h,H);
result=test_close(result,'Hilbert Boltzmann state',rho_h,rho_ref,1e-13,1e-13,...
                  'Hilbert-space thermal equilibrium is exp(-beta H)/Tr[exp(-beta H)]');
result=test_close(result,'Hilbert trace normalisation',trace(rho_h),1,1e-13,1e-13,...
                  'a density matrix returned by equilibrium() must have unit trace');
result=test_close(result,'Hilbert Hermiticity',rho_h,rho_h',1e-13,1e-13,...
                  'a Boltzmann density matrix of a Hermitian Hamiltonian is Hermitian');

% Build the matching Zeeman-Liouville system
bas.formalism='zeeman-liouv';
spin_l=test_spin_system(sys,inter,bas);

% Compare left-product Liouville equilibrium to the vectorised density matrix
H_left=kron(speye(2),H);
rho_l=equilibrium(spin_l,H_left);
result=test_close(result,'Zeeman-Liouville Boltzmann state',rho_l,rho_ref(:),1e-13,1e-13,...
                  'left-product Liouville equilibrium must match the vectorised Hilbert density matrix');

% Build a minimal anisotropic Hamiltonian cell for the orientation branch
Q{1}=cell(3,3);
[Q{1}{:}]=deal(sparse(2,2));
Q{1}{2,2}=sparse(2*pi*diag([2e11 -2e11]));
euler_angles=[0 0 0];
H_oriented=H+orientation(Q,euler_angles);

% Compare the four-argument branch to the explicitly oriented Hamiltonian
rho_orient=equilibrium(spin_h,H,Q,euler_angles);
rho_direct=equilibrium(spin_h,H_oriented);
result=test_close(result,'oriented equilibrium branch',rho_orient,rho_direct,1e-13,1e-13,...
                  'the oriented branch must thermalise H+orientation(Q,euler_angles)');

end


