% Regression for restored point-to-point flux diagnostics. Syntax:
%
%                result=test_kinetics_flux_report()
%
% Outputs: result - flux reporting and unchanged-generator checks
%
% ilya.kuprov@weizmann.ac.il

function result=test_kinetics_flux_report()

% Register the report and generator assertions with the regression suite
result=new_test_result('kernel/cwdm_flux_report','Flux reporting',...
                      'Flux diagnostics preserve directions, rates, and K.');

% Use asymmetric rates so the original from/to convention is observable
sys.magnet=9.4; sys.isotopes={'1H','1H'};
sys.parallel={'local',1};
inter.chem.flux_rate=[0 3;7 0];
flux_types={'intramolecular','intermolecular'};
bas.formalism='sphten-liouv'; bas.approximation={'none'};
for n=1:numel(flux_types)

    % Creation prints the original rate table; kinetics remains quiet
    inter.chem.flux_type=flux_types{n};
    creation=evalc('spin_system=create(sys,inter);');
    spin_system=basis(spin_system,bas);
    ordinary=kinetics(spin_system);
    explicit=evalc('reported=kinetics(spin_system);');
    assert(isempty(explicit),'kinetics must not print summary reports.');
    assert(isequal(ordinary,reported),'kinetics generator changed.');
    text=creation;
    assert(isscalar(strfind(text,'point-to-point flux rates:')),...
           'expected one point-to-point flux table.');
    assert(contains(text,[flux_types{n} ' magnetisation flux']),...
           'magnetisation flux type is missing.');
    assert(~isempty(regexp(text,'1\s+2\s+\+3\.000e\+00','once')),...
           'flux from spin 1 to spin 2 is missing or reversed.');
    assert(~isempty(regexp(text,'2\s+1\s+\+7\.000e\+00','once')),...
           'flux from spin 2 to spin 1 is missing or reversed.');
    fprintf('FLUX_REPORT_CASE %s creation=1 kinetics_quiet=1 generator_unchanged=1\n',flux_types{n});

end

% Absent and all-zero flux retain the original quiet report behavior
for n=1:2
    if n==1
        no_flux=struct();
    else
        no_flux.chem.flux_rate=zeros(2);
        no_flux.chem.flux_type='intramolecular';
    end
    creation=evalc('spin_system=create(sys,no_flux);');
    assert(~contains(creation,'point-to-point flux rates:'));
    spin_system=basis(spin_system,bas);
    ordinary=kinetics(spin_system);
    explicit=evalc('reported=kinetics(spin_system);');
    assert(isempty(explicit),'no-flux report behavior changed.');
    assert(isequal(ordinary,reported),'no-flux generator changed.');
    fprintf('FLUX_REPORT_NO_FLUX case=%d quiet=1 generator_unchanged=1\n',n);
end

disp('TEST_KINETICS_FLUX_REPORT_SUCCESS');

end


