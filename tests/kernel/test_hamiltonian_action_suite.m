% Tests descriptor-backed Hamiltonian action objects. Syntax:
%
%                    result=test_hamiltonian_action_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks that hamiltonian_action objects reproduce the action of
% I+orientation(Q,euler_angles) without assembling that Hamiltonian inside
% the overload path.
%
% ilya.kuprov@weizmann.ac.il
% aditya.dev@weizmann.ac.il

function result=test_hamiltonian_action_suite()

% Announce the test target
fprintf('TESTING: Hamiltonian action overload\n');

% State the Hamiltonian action target of the test
result=new_test_result('kernel/hamiltonian_action_suite',...
                       'Hamiltonian action overload',...
                       'descriptor-backed Hamiltonian action must match explicit oriented Hamiltonians.');

% Check Zeeman Liouville-space action
result=local_case(result,'zeeman-liouv','none',false);

% Check reduced spherical-tensor Liouville-space action
result=local_case(result,'sphten-liouv','IK-0',true);

% Check giant spin terms in Zeeman Liouville space
result=local_giant_case(result);

end




function result=local_case(result,formalism,approximation,reduced)

% Build a compact anisotropic electron-proton system
sys.magnet=1;
sys.isotopes={'E','1H'};
inter.zeeman.matrix={diag([2.00 2.05 2.10]),zeros(3)};
inter.coupling.matrix=cell(2,2);
inter.coupling.matrix{1,2}=2*pi*diag([1.0e6 1.2e6 1.5e6]);
inter.coupling.matrix{2,1}=inter.coupling.matrix{1,2}.';
bas.formalism=formalism;
bas.approximation=approximation;
if reduced
    bas.level=1;
end
spin_system=test_spin_system(sys,inter,bas);
spin_system=assume(spin_system,'labframe');

% Build explicit and descriptor-backed Hamiltonians
euler_angles=[0.31 0.37 0.43];
[I,Q,descr]=hamiltonian(spin_system,'comm');
H_ref=I+orientation(Q,euler_angles);
H_act=hamiltonian_action(spin_system,descr,euler_angles,'comm');

% Compare action on deterministic dense right-hand sides
rng(260602169);
x=randn(size(H_ref,2),3)+1i*randn(size(H_ref,2),3);
test_name=['Hamiltonian action ' formalism '/' approximation];
result=test_close(result,test_name,H_act*x,H_ref*x,1e-9,1e-12,...
                  'descriptor-backed action must match the explicit oriented Hamiltonian');

% Compare sparse conversion on the same compact system
test_name=['Hamiltonian sparse ' formalism '/' approximation];
result=test_close(result,test_name,sparse(H_act),H_ref,1e-9,1e-12,...
                  'sparse conversion must match the explicit oriented Hamiltonian');

% Compare full conversion on the same compact system
test_name=['Hamiltonian full ' formalism '/' approximation];
result=test_close(result,test_name,full(H_act),full(H_ref),1e-9,1e-12,...
                  'full conversion must match the explicit oriented Hamiltonian');

% Check matrix-like predicates and dimensions
if (~isequal(size(H_act),size(H_ref)))||...
   (size(H_act,1)~=size(H_ref,1))||...
   (size(H_act,2)~=size(H_ref,2))||...
   (size(H_act,3)~=1)||...
   (numel(H_act)~=prod(size(H_ref)))||...
   (~ismatrix(H_act))||(~isnumeric(H_act))||iseye(H_act)
    error('FAILED: hamiltonian_action matrix predicates failed.');
end
[dim_a,dim_b,dim_c]=size(H_act);
if (dim_a~=size(H_ref,1))||(dim_b~=size(H_ref,2))||(dim_c~=1)
    error('FAILED: hamiltonian_action size trailing dimensions failed.');
end
result.messages{end+1}=['PASS: hamiltonian_action size and predicates ' formalism '/' approximation '.'];

% Check matrix norm overload
result=test_close(result,['Hamiltonian norm ' formalism '/' approximation],...
                  norm(H_act,1),norm(H_ref,1),1e-9,1e-12,...
                  'norm overload must match the explicit oriented Hamiltonian');

% Check cheap norm dispatch
result=test_close(result,['Hamiltonian cheap norm ' formalism '/' approximation],...
                  cheap_norm(H_act),norm(H_ref,1),1e-9,1e-12,...
                  'cheap_norm must dispatch through the norm overload');

% Check direct descriptor-term operator generation
if height(descr)>0
    xyz=term_oper(H_act,1);
    if (~isnumeric(xyz))||(size(xyz,2)~=3)
        error('FAILED: hamiltonian_action term_oper returned an invalid XYZ array.');
    end
end
result.messages{end+1}=['PASS: hamiltonian_action term operator ' formalism '/' approximation '.'];

end




function result=local_giant_case(result)

% Build a compact giant spin system
sys.magnet=1;
sys.isotopes={'E3'};
inter.zeeman.matrix={diag([2.00 2.01 2.02])};
inter.giant.coeff={{[0 0 0],[0.11e6 -0.03e6 0.07e6 0.02e6 -0.05e6]}};
inter.giant.euler={{[0 0 0],[0.17 0.23 0.31]}};
bas.formalism='zeeman-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);
spin_system=assume(spin_system,'labframe');

% Build explicit and descriptor-backed Hamiltonians
euler_angles=[0.41 0.29 0.13];
[I,Q,descr]=hamiltonian(spin_system,'comm');
H_ref=I+orientation(Q,euler_angles);
H_act=hamiltonian_action(spin_system,descr,euler_angles,'comm');

% Compare action on deterministic dense right-hand sides
rng(260607951);
x=randn(size(H_ref,2),2)+1i*randn(size(H_ref,2),2);
result=test_close(result,'Hamiltonian action giant spin',H_act*x,H_ref*x,...
                  1e-9,1e-12,...
                  'descriptor-backed action must include giant spin terms');

% Compare sparse conversion on the same compact system
result=test_close(result,'Hamiltonian sparse giant spin',sparse(H_act),H_ref,...
                  1e-9,1e-12,...
                  'sparse conversion must include giant spin terms');

end

