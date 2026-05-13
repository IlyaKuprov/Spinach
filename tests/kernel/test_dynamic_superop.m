% Tests superop() sparse-XYZ spherical-tensor product operators. Syntax:
%
%                    result=test_dynamic_superop()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks the unit-operator shortcut, direct commutator and
% anticommutator identities, and Lz spherical-tensor projection eigenvalues.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_superop()

% Announce the test target
fprintf('TESTING: Spherical-tensor superoperator construction\n');

% State the superop target of the test
result=new_test_result('kernel/dynamic_superop',...
                       'Spherical-tensor superoperator construction',...
                       'superop() must produce sparse XYZ product superoperators consistent with angular-momentum algebra.');

% Build a two-spin spherical-tensor Liouville-space system
sys.magnet=0;
sys.isotopes={'1H','13C'};
inter.zeeman.scalar={0,0};
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);
matrix_dim=size(spin_system.bas.basis,1);

% Check the inactive-spin shortcut returns the full identity
A_unit=local_xyz_to_sparse(superop(spin_system,[0 0],'left'),matrix_dim);
result=test_close(result,'superop unit shortcut',A_unit,speye(matrix_dim),1e-15,1e-15,...
                  'an all-zero opspec must map to the unit operator over the full basis');

% Build left, right, commutator, and anticommutator forms for Lz on spin 1
A_left=local_xyz_to_sparse(superop(spin_system,[2 0],'left'),matrix_dim);
A_right=local_xyz_to_sparse(superop(spin_system,[2 0],'right'),matrix_dim);
A_comm=local_xyz_to_sparse(superop(spin_system,[2 0],'comm'),matrix_dim);
A_acomm=local_xyz_to_sparse(superop(spin_system,[2 0],'acomm'),matrix_dim);

% Check algebraic side-product identities
result=test_close(result,'superop comm identity',A_comm,A_left-A_right,1e-15,1e-15,...
                  'a commutator superoperator must equal left multiplication minus right multiplication');
result=test_close(result,'superop acomm identity',A_acomm,A_left+A_right,1e-15,1e-15,...
                  'an anticommutator superoperator must equal left multiplication plus right multiplication');

% Check the Lz commutator eigenvalues on irreducible tensor projections
[~,m_proj]=lin2lm(spin_system.bas.basis(:,1));
result=test_close(result,'superop Lz projection eigenvalues',diag(A_comm),m_proj,1e-15,1e-15,...
                  'the commutator [Lz,T(l,m)] must return m*T(l,m) on the active spin');

% Exercise a two-spin product operator through the multi-spin path
A_left=local_xyz_to_sparse(superop(spin_system,[2 2],'left'),matrix_dim);
A_right=local_xyz_to_sparse(superop(spin_system,[2 2],'right'),matrix_dim);
A_comm=local_xyz_to_sparse(superop(spin_system,[2 2],'comm'),matrix_dim);
result=test_close(result,'superop product comm identity',A_comm,A_left-A_right,1e-15,1e-15,...
                  'multi-spin product commutators must obey the same sided-product identity');
result=test_true(result,'superop product sparsity',nnz(A_comm)>0&&nnz(A_comm)<matrix_dim^2,...
                 'the two-spin product superoperator must be non-zero and sparse in the tensor basis');

end


function A=local_xyz_to_sparse(xyz,matrix_dim)

% Convert Spinach XYZ sparse triples into Matlab sparse storage
A=sparse(xyz(:,1),xyz(:,2),complex(xyz(:,3)),matrix_dim,matrix_dim);

end


