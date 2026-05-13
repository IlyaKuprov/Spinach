% Tests remaining regularisation and inverse-problem utilities. Syntax:
%
%                    result=test_dynamic_remaining_regularisation_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks L-curve corner detection, positivity-constrained
% Tikhonov inversion, and L1 sparsity targeting on compact analytical
% inverse problems.
%
% ilya.kuprov@weizmann.ac.il

function result=test_dynamic_remaining_regularisation_suite()

% State the utility target of the test
result=new_test_result('kernel/dynamic_remaining_regularisation_suite',...
                       'Remaining regularisation utilities',...
                       'Regularisation helpers must recover compact analytical inverse-problem references.');

% Check L-curve analysis on a synthetic corner near lambda equals one
lam=logspace(-3,3,9);
err=sqrt(1+lam.^2);
reg=sqrt(1+lam.^-2);
lam_opt=lcurve(lam,err,reg,'log');
close all;
result=test_true(result,'lcurve synthetic corner',lam_opt>0.1&&lam_opt<10,...
                 'the maximum-curvature point of the symmetric synthetic L-curve should lie near lambda=1');

% Check positivity-constrained Tikhonov inversion against a scalar analytic solution
K=1;
D=1;
KtK=1;
DtD=1;
H=4;
y=2;
lambda=1;
[x_tikh,err_tikh,reg_tikh]=tikhonov(K,D,KtK,DtD,H,y,lambda);
result=test_close(result,'tikhonov scalar solution',x_tikh,1,1e-8,1e-8,...
                  'minimising (x-2)^2+x^2 with x>=0 gives x=1');
result=test_close(result,'tikhonov scalar error',err_tikh,1,1e-8,1e-8,...
                  'the residual error at x=1 is (1-2)^2');
result=test_close(result,'tikhonov scalar regularisation',reg_tikh,1,1e-8,1e-8,...
                  'the regularisation signal at x=1 is x^2');

% Check L1 sparsity targeting on an identity sensing matrix
A=eye(3);
y=[2;0;0];
[x_l1,err_l1,reg_l1]=tikhol1n(A,y,1);
result=test_true(result,'tikhol1n sparse support',nnz(abs(x_l1)>1e-8)==1&&x_l1(1)>0,...
                 'identity data with one non-zero target should keep only the populated component');
result=test_true(result,'tikhol1n finite metrics',isfinite(err_l1)&&isfinite(reg_l1)&&...
                 err_l1>=0&&reg_l1>0,...
                 'tikhol1n() must return finite non-negative error and positive L1 regularisation metrics');

end


