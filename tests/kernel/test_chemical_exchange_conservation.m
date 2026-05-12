% Tests conservation in two-site chemical exchange. Syntax:
%
%                    result=test_chemical_exchange_conservation()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test builds a symmetric two-site exchange model and checks that the
% kinetics generator conserves the total population over the two sites.
%
% ilya.kuprov@weizmann.ac.il

function result=test_chemical_exchange_conservation()

% State the kinetics target of the test
result=new_test_result('kernel/chemical_exchange_conservation',...
                       'Chemical exchange conservation',...
                       'closed two-site exchange must conserve total spin population.');

% Build a symmetric two-site exchange system
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0 0};
inter.chem.parts={1,2};
inter.chem.rates=[-3 3;3 -3];
inter.chem.concs=[1 1];
bas.formalism='sphten-liouv';
bas.approximation='none';
spin_system=test_spin_system(sys,inter,bas);

% Build the kinetics generator
K=kinetics(spin_system);

% Closed Markov kinetics conserve total population by zero column sums
col_sums=sum(full(K),1);
result=test_close(result,'zero column sums',col_sums,zeros(size(col_sums)),1e-15,1e-15,...
                  'probability leaving each site must enter the other site');

end

