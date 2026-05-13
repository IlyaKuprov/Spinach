% Tests sparse, tensor, and numerical utility helpers. Syntax:
%
%                    result=test_sparse_tensor_utility_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks Kronecker-product application, Blicharski invariants,
% spectral densities, SVD truncation helpers, sparse density, and related
% numerical utilities against direct definitions.
%
% ilya.kuprov@weizmann.ac.il

function result=test_sparse_tensor_utility_suite()

% Announce the test target
fprintf('TESTING: Sparse, tensor, and numerical utility functions\n');

% State the utility target of the test
result=new_test_result('kernel/sparse_tensor_utility_suite',...
                       'Sparse, tensor, and numerical utility functions',...
                       'sparse/tensor utilities must match direct dense algebra on small reference cases.');

% Kronecker-product application must match explicit kron products
Q={sparse([1 2;3 4]),sparse([0 5;6 7])};
x=(1:4)';
K=kron(Q{1},Q{2});
result=test_close(result,'kronm vector product',kronm(Q,x),K*x,1e-14,1e-14,...
                  'kronm applies a Kronecker product without explicitly forming it');
result=test_close(result,'kronm_new vector product',kronm_new(Q,x),K*x,1e-14,1e-14,...
                  'kronm_new applies the same Kronecker product using tensor contractions');

% Blicharski invariants ignore isotropic trace and follow their definitions
A=[1 2 3;4 5 6;7 8 9];
[Lsq,Dsq]=blinv(A);
Lsq_ref=(A(1,2)-A(2,1))^2+(A(1,3)-A(3,1))^2+(A(2,3)-A(3,2))^2;
Dsq_ref=A(1,1)^2+A(2,2)^2+A(3,3)^2-A(1,1)*A(2,2)-A(1,1)*A(3,3)-A(2,2)*A(3,3)+...
        (3/4)*((A(1,2)+A(2,1))^2+(A(1,3)+A(3,1))^2+(A(2,3)+A(3,2))^2);
result=test_close(result,'blinv first-rank invariant',Lsq,Lsq_ref,1e-14,1e-14,...
                  'first-rank Blicharski invariant is the squared antisymmetric tensor amplitude');
result=test_close(result,'blinv second-rank invariant',Dsq,Dsq_ref,1e-14,1e-14,...
                  'second-rank Blicharski invariant is the traceless symmetric tensor amplitude');
B=[2 -1 0; 3 4 1; 5 6 -2];
[X1,X2]=blprod(A,B);
[Lam,Dam]=blinv(A-B); [Lap,Dap]=blinv(A+B);
result=test_close(result,'blprod first-rank polarisation',X1,(Lap-Lam)/4,1e-14,1e-14,...
                  'first-rank tensor product uses the polarisation identity');
result=test_close(result,'blprod second-rank polarisation',X2,(Dap-Dam)/4,1e-14,1e-14,...
                  'second-rank tensor product uses the polarisation identity');

% Lorentzian spectral density follows its closed definition
L=2; Drot=1.5e6; omega=2.0e5; tau=1/(L*(L+1)*Drot);
result=test_close(result,'spden Lorentzian',spden(L,Drot,omega),(tau/(2*L+1))/(1+(tau*omega)^2),1e-20,1e-14,...
                  'spectral density is the rank-dependent Lorentzian rotational correlation function');

% Frobenius SVD truncation keeps the smallest rank whose dropped tail is below tolerance
s=[5 1 0.01];
result=test_close(result,'frob_chop keep count',frob_chop(s,0.02),2,0,0,...
                  'with tolerance 0.02, only the 0.01 singular value may be dropped in Frobenius norm');

% svd_shrink factorisation must reconstruct retained singular content
spin_system.sys.output='hush';
rho=diag([4 1 1e-6]);
[vec,cov]=svd_shrink(spin_system,rho,1e-4);
result=test_close(result,'svd_shrink reconstruction',vec*cov',diag([4 1 0]),1e-12,1e-12,...
                  'svd_shrink returns vector/covector factors for singular values above threshold');

end
