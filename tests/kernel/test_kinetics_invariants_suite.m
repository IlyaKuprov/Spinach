% Tests deterministic chemical kinetics helpers. Syntax:
%
%                    result=test_kinetics_invariants_suite()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test checks closed-form steady states, independent reaction blocks,
% first-order exchange routing and conservation, and empty reaction maps
% on a reordered local descriptor.
%
% ilya.kuprov@weizmann.ac.il

function result=test_kinetics_invariants_suite()

% Announce the test target
fprintf('TESTING: Chemical kinetics invariants\n');

% State the kinetics target of the test
result=new_test_result('kernel/kinetics_invariants_suite',...
                       'Chemical kinetics invariants',...
                       'supported kinetic generators must conserve and route spin order.');

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

% Build a two-substance one-way exchange system
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0,0};
inter.chem.parts={1,2}; inter.chem.concs=[1 1];
inter.chem.reactions={struct('reactants',1,'products',2,...
    'matching',[1 2],'rate',3)};
bas.formalism='sphten-liouv';
bas.approximation={'none','none'};
spin_system=test_spin_system(sys,inter,bas);

% First-order exchange must conserve every column sum of the kinetic generator
K=kinetics(spin_system);
col_sums=sum(full(K),1);
result=test_close(result,'kinetics exchange column sums',col_sums,zeros(size(col_sums)),1e-14,1e-14,...
                  'first-order exchange moves spin order without losing its column sum');

% One-way flux drains source spin order and fills matched destination spin order
rho_source=state(spin_system,'Lz',1);
rho_destin=state(spin_system,'Lz',2);
result=test_close(result,'kinetics source-to-destination routing',K*rho_source,3*(rho_destin-rho_source),1e-14,1e-14,...
                  'flux transfers longitudinal spin order from spin one to spin two at the specified rate');
result=test_close(result,'kinetics no action on destination source',K*rho_destin,zeros(size(rho_destin)),1e-14,1e-14,...
                  'one-way flux does not drain spin order already on its destination');

% Empty reactions traverse a single local descriptor with permuted spin columns
inter.chem.parts={[2 1]}; inter.chem.concs=1; inter.chem.reactions={};
bas.approximation={'none'};
spin_system=test_spin_system(sys,inter,bas);
reaction.reactants=[]; reaction.products=[]; reaction.matching=zeros(0,2);
G=react_gen(spin_system,reaction);
result=test_true(result,'react_gen reordered local columns',...
                 iscell(G)&&isempty(G),...
                 'an empty reaction is valid independently of the local spin-column order');

end


