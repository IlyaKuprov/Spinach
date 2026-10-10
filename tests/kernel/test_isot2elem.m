% Tests isotope-number stripping and cell-array shape preservation.
% Syntax:
%
%                       result=test_isot2elem()
%
% Outputs:
%
%    result  - regression test result with explanatory messages
%
% talos@spindynamics.org

function result=test_isot2elem()

% Announce the test target
fprintf('TESTING: Isotope-to-element conversion\n');
result=new_test_result('kernel/isot2elem','Isotope-to-element conversion',...
                       'Strip isotope numbers without changing cell-array shape or order.');

% Check ordinary nuclei and two-letter element symbols
result=test_true(result,'nuclear labels',...
                 isequal(isot2elem({'1H','13C','15N','35Cl','129Xe'}),...
                                  {'H','C','N','Cl','Xe'}),...
                 'mass numbers must be removed, and element capitalisation retained');

% Check cell-array shapes and repeated labels
result=test_true(result,'column cells',...
                 isequal(isot2elem({'2H';'1H';'13C'}),{'H';'H';'C'}),...
                 'conversion must preserve the column shape and label order');
result=test_true(result,'matrix cells',...
                 isequal(isot2elem({'1H','13C';'35Cl','129Xe'}),...
                                  {'H','C';'Cl','Xe'}),...
                 'conversion must retain both matrix dimensions');
result=test_true(result,'empty cells',...
                 isequal(isot2elem(cell(0,3)),cell(0,3)),...
                 'an empty input must retain its cell-array dimensions');

% Check that conversion is literal rather than an isotope lookup
result=test_true(result,'non-digit characters',...
                 isequal(isot2elem({'H','99Tc_m'}),{'H','Tc_m'}),...
                 'the function removes digits only, not isotope-state suffixes');

% Check rejection of unsupported input types and character matrices
bad_inputs={'1H',{'1H',13},{['1H';'2H']},{{'1H'}}};
for n=1:numel(bad_inputs)
    caught=false();
    try
        isot2elem(bad_inputs{n});
    catch fault
        caught=strcmp(fault.message,...
                      'isotopes must be a cell array of character row vectors.');
    end
    result=test_true(result,['invalid input ' num2str(n)],caught,...
                     'unsupported inputs must receive the documented type diagnostic');
end

end


