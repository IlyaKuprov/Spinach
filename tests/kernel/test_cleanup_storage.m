% Tests storage policy after fused matrix multiplication. Syntax:
%
%                    result=test_cleanup_storage()
%
% Outputs:
%
%     result - regression result with rounding and representation checks
%
% ilya.kuprov@weizmann.ac.il

function result=test_cleanup_storage()

% State the test target and ordinary cleanup policy
result=new_test_result('kernel/cleanup_storage','Product cleanup policy',...
                       'Cleanup controls final storage; roundoff is permitted.');
s.sys.disable={}; s.sys.parallel={'processes',1}; s.tols.small_matrix=10; s.tols.dense_matrix=0.5;

% Check sparse and dense outcomes without changing the supplied tolerances
for density=[0.1 0.9]
    A=sprand(100,100,density); B=speye(100); tol=1e-6;
    reference=clean_up(s,A*B,tol); observed=multiclean(s,A,B,tol);
    result=test_true(result,'numerical agreement',...
                     norm(observed-reference,'fro')<=...
                     100*eps*max(1,norm(reference,'fro')),...
                     'Panel cleanup must agree within numerical roundoff.');
    result=test_true(result,'global storage',issparse(observed)==issparse(reference),...
                     'Final storage must follow clean_up.');
end

% Zero and disabled cleanup are handled by clean_up itself
A=randn(12); B=randn(12);
for disabled=[false true]
    if disabled, s.sys.disable={'clean-up'}; end
    observed=multiclean(s,A,B,0); reference=A*B;
    result=test_true(result,'unrounded product',...
                     norm(observed-reference,'fro')<=12*eps*norm(reference,'fro'),...
                     'Unrounded products must agree within numerical roundoff.');
end

end


