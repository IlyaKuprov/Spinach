% Tests removal of the isotropic tensor trace. Syntax:
%
%                    result=test_remtrace_tensor()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that remtrace subtracts the isotropic component of a
% second-rank interaction tensor, leaving the anisotropic traceless part.
%
% ilya.kuprov@weizmann.ac.il

function result=test_remtrace_tensor()

% Announce the test target
fprintf('TESTING: Traceless rank-two tensor construction\n');

% State the physical target of the test
result=new_test_result('kernel/remtrace_tensor',...
                       'Traceless rank-two tensor construction',...
                       'anisotropic interaction tensors are obtained by subtracting the isotropic trace.');

% Define a symmetric interaction tensor with non-zero isotropic part
A=[1 2 0;2 3 0;0 0 5];
A_ref=A-eye(3)*trace(A)/3;
A_obs=remtrace(A);

% Check explicit trace removal and invariants
result=test_close(result,'explicit isotropic subtraction',A_obs,A_ref,1e-15,1e-15,...
                  'the isotropic part is trace(A)/3 times the unit matrix');
result=test_close(result,'zero trace',trace(A_obs),0,1e-15,1e-15,...
                  'anisotropic second-rank tensors have zero trace');
result=test_close(result,'anisotropy preserved',A_obs(1,2),A(1,2),1e-15,1e-15,...
                  'subtracting isotropic trace leaves off-diagonal anisotropy unchanged');

end

