% Tests linear-algebra, angular-momentum, and perturbation utilities. Syntax:
%
%                    result=test_linear_perturbation_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks spin-addition projectors, Rayleigh-Schrödinger
% and Van Vleck perturbation theory, analytical Tikhonov inversion,
% transfer matrices, and finite-difference Jacobian estimation.
%
% ilya.kuprov@weizmann.ac.il

function result=test_linear_perturbation_suite()

% Announce the test target
fprintf('TESTING: Linear algebra and perturbation utilities\n');

% State the utility target of the test
result=new_test_result('kernel/linear_perturbation_suite',...
                       'Linear algebra and perturbation utilities',...
                       'Small analytical linear-algebra cases must match exact matrix and perturbation references.');

% Check spin-half addition into singlet and triplet irreducible blocks
[multiplicities,projectors]=add_spins(1/2,1/2);
result=test_close(result,'add_spins multiplicities',multiplicities,[1 3],...
                  1e-14,1e-14,...
                  'two spin-half irreps reduce into one singlet and one triplet');
projector_sum=projectors{1}*projectors{1}'+projectors{2}*projectors{2}';
result=test_close(result,'add_spins projector completeness',projector_sum,eye(4),...
                  1e-13,1e-13,...
                  'singlet and triplet projectors resolve the full product space');
result=test_close(result,'add_spins singlet orthonormality',...
                  projectors{1}'*projectors{1},eye(1),1e-13,1e-13,...
                  'the singlet projector columns are orthonormal');
result=test_close(result,'add_spins triplet orthonormality',...
                  projectors{2}'*projectors{2},eye(3),1e-13,1e-13,...
                  'the triplet projector columns are orthonormal');

% Check second-order perturbation energy shifts for a two-level system
base_energy=[0;10];
pert_strength=0.01;
pert_mat=[0 pert_strength;pert_strength 0];
pert_ref=[-pert_strength^2/10;10+pert_strength^2/10];
[rs_energy,rs_vectors]=rspert(base_energy,pert_mat,2);
result=test_close(result,'rspert two-level energy shifts',rs_energy,pert_ref,...
                  1e-13,1e-13,...
                  'second-order RSPT gives the textbook avoided-crossing energy shifts');
result=test_close(result,'rspert vector normalisation',...
                  sum(abs(rs_vectors).^2,1),ones(1,2),1e-14,1e-14,...
                  'RSPT eigenvectors are column-normalised');

% Check Van Vleck perturbation theory on the same two-level system
[vv_energy,vv_gen]=vvpert(base_energy,pert_mat,2);
result=test_close(result,'vvpert two-level energy shifts',vv_energy,pert_ref,...
                  1e-13,1e-13,...
                  'second-order Van Vleck perturbation theory gives the same energy shifts');
result=test_close(result,'vvpert anti-Hermitian generator',vv_gen+vv_gen',zeros(2),...
                  1e-14,1e-14,...
                  'the real Van Vleck generator is antisymmetric for a Hermitian perturbation');

% Check the analytical indefinite Tikhonov solution for identity operators
fit_mat=eye(2);
reg_mat=eye(2);
fit_rhs=[3;6];
reg_param=1/2;
[tikh_x,tikh_err,tikh_reg]=tikhoind(fit_mat,reg_mat,fit_rhs,reg_param);
tikh_ref=fit_rhs/(1+reg_param);
result=test_close(result,'tikhoind identity solution',tikh_x,tikh_ref,...
                  1e-14,1e-14,...
                  'identity data and regularisation matrices shrink the right-hand side by 1+lambda');
result=test_close(result,'tikhoind error signal',tikh_err,norm(tikh_ref-fit_rhs,2)^2,...
                  1e-14,1e-14,...
                  'the reported fit error is norm(K*x-y)^2');
result=test_close(result,'tikhoind regularisation signal',tikh_reg,norm(tikh_ref,2)^2,...
                  1e-14,1e-14,...
                  'the reported regularisation signal is norm(D*x)^2');

% Check transfer matrix recovery from more vector pairs than dimensions
amp_inputs=[1 0 1 2;0 1 1 -1];
transfer_ref=[2 -1;1/2 3];
amp_outputs=transfer_ref*amp_inputs;
transfer_obs=transfermat(amp_inputs,amp_outputs);
result=test_close(result,'transfermat exact recovery',transfer_obs,transfer_ref,...
                  1e-13,1e-13,...
                  'a full-row-rank input stack recovers the exact linear transfer matrix');

% Check finite-difference Jacobian estimation against analytical derivatives
jac_fun=@(x)[x(1)^2+3*x(2);sin(x(1)*x(2))];
jac_point=[2;0.3];
[jac_obs,jac_err]=jacobianest(jac_fun,jac_point);
jac_ref=[4 3;0.3*cos(0.6) 2*cos(0.6)];
result=test_close(result,'jacobianest analytical Jacobian',jac_obs,jac_ref,...
                  1e-6,1e-6,...
                  'the finite-difference Jacobian matches the analytical derivative matrix');
result=test_true(result,'jacobianest finite error estimates',...
                 all(isfinite(jac_err(:)))&&all(jac_err(:)>=0),...
                 'Jacobian error estimates are finite and non-negative');

end


