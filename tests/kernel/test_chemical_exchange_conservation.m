% Tests two-site chemical exchange conservation. Syntax:
%
%                    result=test_chemical_exchange_conservation()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The symmetric two-site fixture checks matter and longitudinal spin
% conservation, matched transfer, and the analytic population trajectory.
%
% ilya.kuprov@weizmann.ac.il

function result=test_chemical_exchange_conservation()

% Announce the test target
fprintf('TESTING: Chemical exchange conservation\n');

% State the kinetics target of the test
result=new_test_result('kernel/chemical_exchange_conservation',...
                       'Chemical exchange conservation',...
                       'two-site reaction records must conserve matter and transfer matched spin order.');

% Build a symmetric two-site exchange system
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0 0};
inter.chem.parts={1,2};
inter.chem.reactions={struct('reactants',1,'products',2,...
    'matching',[1 2],'rate',3),...
    struct('reactants',2,'products',1,...
    'matching',[2 1],'rate',3)};
inter.chem.concs=[1 1];
bas.formalism='sphten-liouv';
bas.approximation={'none','none'};
spin_system=test_spin_system(sys,inter,bas);

% Check total concentration and longitudinal spin conservation
K=kinetics(spin_system);
units=spin_system.bas.offsets(1:end-1)+1;
result=test_close(result,'concentration conservation',sum(K(units,:),1),...
                  zeros(1,size(K,2)),1e-14,1e-14,'total molecular population is conserved');
coil=coil_state(spin_system,'Lz','1H','exact');
result=test_close(result,'longitudinal conservation',coil'*K,...
                  zeros(1,size(K,2)),1e-14,1e-14,'matched longitudinal spin order is conserved');

% Check matched spin transfer and the analytic population trajectory
rho_a=state(spin_system,'Lz',1);
rho_b=state(spin_system,'Lz',2);
result=test_close(result,'matched transfer',K*rho_a,3*(rho_b-rho_a),...
                  1e-14,1e-14,'the first-order record transfers spin order at rate three');
rho=unit_state(spin_system); rho(units)=[2;0];
rho=expm(full(K)*0.2)*rho;
result=test_close(result,'population trajectory',rho(units),...
                  expm([-3 3;3 -3]*0.2)*[2;0],1e-13,1e-13,...
                  'unit coordinates follow the independent two-site rate matrix');

end

