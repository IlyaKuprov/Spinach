% Tests the commutator utility. Syntax:
%
%                    result=test_commutator_utility()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that comm(A,B) implements AB-BA and returns zero for
% mutually commuting matrices.
%
% ilya.kuprov@weizmann.ac.il

function result=test_commutator_utility()

% Announce the test target
fprintf('TESTING: Matrix commutator utility\n');

% State the mathematical target of the test
result=new_test_result('kernel/commutator_utility',...
                       'Matrix commutator utility',...
                       'comm(A,B) must return AB-BA exactly.');

% Define a non-commuting pair and its reference commutator
A=[1 2;3 4];
B=[0 1;-1 2];
C=A*B-B*A;

% Check the utility on a non-commuting pair
result=test_close(result,'non-commuting pair',comm(A,B),C,1e-15,1e-15,...
                  'the commutator is defined as AB-BA');

% Check a commuting diagonal pair
D=diag([1 2 3]);
E=diag([4 5 6]);
result=test_close(result,'commuting diagonal pair',comm(D,E),zeros(3),1e-15,1e-15,...
                  'diagonal matrices in the same basis commute');

end

