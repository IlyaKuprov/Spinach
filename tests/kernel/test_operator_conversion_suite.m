% Tests Hilbert-to-Liouville operator conversion utilities. Syntax:
%
%                    result=test_operator_conversion_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks direct vectorisation identities for left, right,
% commutation, and anticommutation superoperators, unit_oper dimensions,
% and Lindbladian rate calibration.
%
% ilya.kuprov@weizmann.ac.il

function result=test_operator_conversion_suite()

% Announce the test target
fprintf('TESTING: Hilbert/Liouville conversion helpers\n');

% State the operator-conversion target of the test
result=new_test_result('kernel/operator_conversion_suite',...
                       'Hilbert/Liouville conversion helpers',...
                       'operator conversion functions must implement vectorised product identities.');

% Define a non-trivial Hermitian operator
S=pauli(2);
H=S.z+0.2*S.x;
unit=speye(2);

% Check direct vectorisation formulas
result=test_close(result,'hilb2liouv left',hilb2liouv(H,'left'),kron(unit,H),1e-15,1e-15,...
                  'left multiplication vectorises as kron(I,H)');
result=test_close(result,'hilb2liouv right',hilb2liouv(H,'right'),kron(transpose(H),unit),1e-15,1e-15,...
                  'right multiplication vectorises as kron(transpose(H),I)');
result=test_close(result,'hilb2liouv comm',hilb2liouv(H,'comm'),kron(unit,H)-kron(transpose(H),unit),1e-15,1e-15,...
                  'a commutator is left multiplication minus right multiplication');
result=test_close(result,'hilb2liouv acomm',hilb2liouv(H,'acomm'),kron(unit,H)+kron(transpose(H),unit),1e-15,1e-15,...
                  'an anticommutator is left multiplication plus right multiplication');
result=test_close(result,'hilb2liouv statevec',hilb2liouv(H,'statevec'),H(:),1e-15,1e-15,...
                  'state-vector conversion stacks matrix columns');

% Check unit operator dimensions in major formalisms
sys.magnet=0;
sys.isotopes={'1H'};
inter.zeeman.scalar={0};
bas.formalism='zeeman-hilb';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);
result=test_close(result,'zeeman-hilb unit_oper',unit_oper(spin_system),speye(2),1e-15,1e-15,...
                  'Hilbert-space unit dimension is the spin multiplicity product');
bas.formalism='zeeman-liouv';
spin_system=test_spin_system(sys,inter,bas);
result=test_close(result,'zeeman-liouv unit_oper',unit_oper(spin_system),speye(4),1e-15,1e-15,...
                  'Zeeman Liouville unit dimension is the square of Hilbert dimension');

% Check Lindbladian calibration on a simple vector
A_left=diag([1 0]);
A_right=diag([0 1]);
rho=[1;1];
rate=3.5;
R=lindbladian(A_left,A_right,rho,rate);
obs=real((rho'*R*rho)/(rho'*rho));
result=test_close(result,'lindbladian rate calibration',obs,-rate,1e-12,1e-12,...
                  'lindbladian() rescales the generator to the requested experimental decay rate');

end

