% Tests deterministic chemical kinetics helpers. Syntax:
%
%                    result=test_kinetics_invariants_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks closed-form steady states, independent reaction blocks,
% reaction-generator state routing, and conservation in a tiny exchange
% kinetics superoperator.
%
% ilya.kuprov@weizmann.ac.il

function result=test_kinetics_invariants_suite()

% State the kinetics target of the test
result=new_test_result('kernel/kinetics_invariants_suite',...
                       'Chemical kinetics invariants',...
                       'kinetic generators must conserve matter and route spin order between declared species.');

% Check a two-site steady state from detailed balance
kf=2;
kr=5;
K=[-kf kr; kf -kr];
c0=[2;1];
ctot=sum(c0);
c_ref=ctot*[kr; kf]/(kf+kr);
result=test_close(result,'equilibrate two-site detailed balance',equilibrate(K,c0),c_ref,1e-13,1e-13,...
                  'at equilibrium k_forward c_1 equals k_reverse c_2 and total concentration is conserved');

% Check recursive treatment of independent reaction blocks
K1=[-1 4;1 -4];
K2=[-3 2;3 -2];
K=blkdiag(K1,K2);
c0=[3;0;1;2];
c_ref=[sum(c0(1:2))*[4;1]/5; sum(c0(3:4))*[2;3]/5];
result=test_close(result,'equilibrate independent blocks',equilibrate(K,c0),c_ref,1e-13,1e-13,...
                  'independent kinetic components equilibrate separately and retain their own material totals');

% Check the zero-concentration shortcut
result=test_close(result,'equilibrate zero concentration',equilibrate(K,zeros(4,1)),zeros(4,1),1e-15,1e-15,...
                  'a zero initial concentration vector remains zero for linear kinetics');

% Build a two-site spherical-tensor spin system for exchange and reactions
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0,0};
inter.chem.parts={1,2};
inter.chem.rates=[-3 3;3 -3];
inter.chem.concs=[1 1];
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Closed exchange must conserve every column sum of the kinetic generator
K=kinetics(spin_system);
col_sums=sum(full(K),1);
result=test_close(result,'kinetics closed-exchange column sums',col_sums,zeros(size(col_sums)),1e-14,1e-14,...
                  'in a closed two-site exchange model all probability leaving a column re-enters elsewhere');

% A one-way reaction generator must drain source spin order and fill matched product spin order
reaction.reactants=1;
reaction.products=2;
reaction.matching=[1 2];
G=react_gen(spin_system,reaction);
rho_source=state(spin_system,'Lz',1);
rho_destin=state(spin_system,'Lz',2);
result=test_close(result,'react_gen source-to-product routing',G{1}*rho_source,rho_destin-rho_source,1e-14,1e-14,...
                  'a matched reactant spin order is removed from the source species and inserted into the product species');
result=test_close(result,'react_gen no action on product source',G{1}*rho_destin,zeros(size(rho_destin)),1e-14,1e-14,...
                  'the reactant generator does not drain states already located on the product species');

end


