% Tests global storage conversion after fused matrix multiplication. Syntax:
%
%                    result=test_cleanup_storage()
%
% Outputs:
%
%     result - regression result with rounding and representation checks
%
% ilya.kuprov@weizmann.ac.il

function result=test_cleanup_storage()

% State the test target
result=new_test_result('kernel/cleanup_storage','Exactly-once product cleanup',...
                       'Storage conversion must not round a fused product twice.');

% Distinguish one rounding operation from two at large dynamic range
s.sys.disable={}; s.tols.small_matrix=0; s.tols.dense_matrix=0;
x=-302204600072882.88; tol=0.1;
once=clean_up(s,x,tol); twice=clean_up(s,once,tol);
result=test_true(result,'discriminating rounding case',~isequal(once,twice),...
                 'This input must expose the non-idempotence of floating-point rounding.');
result=test_true(result,'round once',isequal(multiclean(s,1,x,tol),once),...
                 'Fused multiplication must round exactly once.');

% Global sparse conversion must not alter the once-rounded values
s.tols.dense_matrix=1;
A=[x 0;0 0]; expected=clean_up(s,A,tol);
C=multiclean(s,eye(2),A,tol);
result=test_true(result,'sparse conversion',issparse(C)&&isequal(C,expected),...
                 'Global sparse conversion must preserve once-rounded values.');
s.tols.small_matrix=3;
C=multiclean(s,speye(2),sparse(A),tol);
result=test_true(result,'full conversion',~issparse(C)&&isequal(C,full(expected)),...
                 'Small nonzero products must become full without another rounding.');

% The disabled and NaN branches preserve native values and representation
s.sys.disable={'clean-up'};
result=test_true(result,'disabled cleanup',isequal(multiclean(s,1,x,tol),x),...
                 'Disabled cleanup must return the native product.');
s.sys.disable={};
result=test_true(result,'NaN cleanup',isequal(multiclean(s,1,x,NaN),x),...
                 'NaN tolerance must return the native product.');

end


