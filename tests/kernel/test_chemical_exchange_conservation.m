% Tests the two-site chemical exchange implementation boundary. Syntax:
%
%                    result=test_chemical_exchange_conservation()
%
% Outputs:
%
%     result  - regression test result with explanatory messages
%
% The test retains a symmetric two-site exchange fixture and checks the
% named rejection until segmented reaction-record chemistry is available.
% It does not assert numerical conservation for unsupported chemistry.
%
% ilya.kuprov@weizmann.ac.il

function result=test_chemical_exchange_conservation()

% Announce the test target
fprintf('TESTING: Chemical exchange implementation boundary\n');

% State the kinetics target of the test
result=new_test_result('kernel/chemical_exchange_conservation',...
                       'Chemical exchange implementation boundary',...
                       'unsupported two-site exchange must raise the named chemistry boundary error.');

% Build a symmetric two-site exchange system
sys.magnet=14.1;
sys.isotopes={'1H','1H'};
inter.zeeman.scalar={0 0};
inter.chem.parts={1,2};
inter.chem.rates=[-3 3;3 -3];
inter.chem.concs=[1 1];
bas.formalism='sphten-liouv';
bas.approximation={'none','none'};
spin_system=test_spin_system(sys,inter,bas);

% Require the explicit boundary rather than a numerical exchange generator
rejected=false;
try
    kinetics(spin_system);
catch err
    rejected=strcmp(err.identifier,'Spinach:kinetics:segmentedChemistry')&&...
             contains(err.message,'reaction-record implementation (WP3)');
end
result=test_true(result,'kinetics WP3 boundary',rejected,...
                 'unsupported multi-substance chemistry raises the named boundary error');

end

